# -*- coding: utf-8 -*-
"""
FULL_CATALOG_SCAN_PLOT.py (Updated with Bayesian Evidence & Bayes Factor Calculations)

Combined script:
  - CSV-based full catalog RMS & Bayes Factor pipeline (GR vs PQG-PN)
  - Per-event plots (time domain + RMS(beta2) + Bayes factor evaluation)
  - Catalog-level plots:
        * beta2 histograms per approximant
        * Δ distribution (RMS improvement) + normal fit
        * ln B (log-Bayes factor) distribution across catalog
        * cumulative Δ(N) and cumulative ln B(N)
        * Z-score growth Z(N)

Outputs:
  - summary.txt (per-event results including log_L, log_Z, and ln_B)
  - summary_clean.txt (NaN-free global statistics)
  - plots_catalog/<event>_<approximant>_GR_vs_PQGPN.png
  - plots_catalog/<event>_<approximant>_RMS_beta2.png
  - plots_catalog/beta2_posterior_<approximant>.png
  - plots_catalog/delta_distribution.png
  - plots_catalog/bayes_factor_distribution.png
  - plots_catalog/cumulative_delta.png
  - plots_catalog/zscore_growth.png
"""

import os
import csv
from datetime import datetime
from math import erf, sqrt

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from gwpy.timeseries import TimeSeries
from pycbc.waveform import get_td_waveform
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count

print("[BOOT] Module import completed (including scipy.integrate.simpson for Bayesian evidence)")

# =========================================================
# 0) PARAMETERS
# =========================================================

APPROXIMANTS = [
    "IMRPhenomD",
    "SEOBNRv4",
]

print(f"[CONFIG] Approximants: {APPROXIMANTS}")

detectors = ["H1", "L1"]
print(f"[CONFIG] Detectors: {detectors}")

whiten_seg = 4.0
f_low = 30.0
f_high = 300.0

fit_tmin = -0.10
fit_tmax =  0.02

beta2_min = -1.0e-4
beta2_max =  1.0e-4
beta2_N   = 41

N_b1 = 21
N_b0 = 21

b1_min = -5e-3
b1_max =  5e-3

b0_min = -np.pi
b0_max =  np.pi

print(f"[CONFIG] beta2 grid: [{beta2_min}, {beta2_max}] with N={beta2_N}")
print(f"[CONFIG] b1 grid: [{b1_min}, {b1_max}] with N={N_b1}")
print(f"[CONFIG] b0 grid: [{b0_min}, {b0_max}] with N={N_b0}")

# fetch_open_data cache (detector, t_start, t_end) -> TimeSeries
_FETCH_CACHE = {}

# plot output directory
PLOT_DIR = "plots_catalog"
os.makedirs(PLOT_DIR, exist_ok=True)
print(f"[INIT] Plot directory ensured: {PLOT_DIR}")


# =========================================================
# 1) CSV LOADING
# =========================================================

def load_events_from_csv(csv_path="events.csv", max_events=None):
    print(f"[CSV] Loading events from {csv_path} (max_events={max_events})")
    EVENTS = {}
    count = 0

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                name = row["name"]
                gps = float(row["gps"])
                m1 = float(row["mass_1_source"])
                m2 = float(row["mass_2_source"])
            except Exception as e:
                print(f"[CSV][WARN] Skipping malformed row: {row} | error={e}")
                continue

            f_lower = 30.0
            fetch_win = 4.0

            EVENTS[name] = (gps, (m1, m2), f_lower, fetch_win)
            count += 1
            print(f"[CSV] Loaded event {name}: gps={gps}, masses=({m1},{m2})")

            if max_events is not None and count >= max_events:
                print("[CSV] Reached max_events limit")
                break

    print(f"[CSV] Total events loaded: {len(EVENTS)}")
    return EVENTS


# =========================================================
# 2) HELPER & BAYESIAN FUNCTIONS
# =========================================================

def rms_error(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))


def compute_log_likelihood(y_data, h_model):
    """
    Gaussian log-likelihood for whitened strain data with unit noise variance.
    log L = -0.5 * sum((y_data - h_model)^2)
    """
    return -0.5 * np.sum((y_data - h_model)**2)


def fetch_open_data_cached(det, t_start, t_end):
    key = (det, float(t_start), float(t_end))
    if key in _FETCH_CACHE:
        print(f"[CACHE] Hit for {det} {t_start}..{t_end}")
        return _FETCH_CACHE[key]
    print(f"[CACHE] Miss for {det} {t_start}..{t_end} — fetching")
    h = TimeSeries.fetch_open_data(det, t_start, t_end, cache=True)
    _FETCH_CACHE[key] = h
    return h


def load_and_condition(det, t0, fetch_win):
    t_start = t0 - fetch_win/2.0
    t_end   = t0 + fetch_win/2.0

    h = fetch_open_data_cached(det, t_start, t_end)
    white = h.whiten(whiten_seg)
    bp = white.bandpass(f_low, f_high)
    window = bp.crop(t0 - 0.20, t0 + 0.05)

    t_data = window.times.value - t0
    y_data = window.value
    dt = window.dt.value

    mask_fit = (t_data > fit_tmin) & (t_data < fit_tmax)

    return {
        "t_data": t_data,
        "y_data": y_data,
        "dt": dt,
        "mask_fit": mask_fit,
        "t_fit": t_data[mask_fit],
        "y_fit": y_data[mask_fit],
    }


def generate_gr_on_data_grid(approximant, dt_data, t_data, m1, m2, f_lower):
    hp, hc = get_td_waveform(
        approximant=approximant,
        mass1=m1,
        mass2=m2,
        spin1z=0.0,
        spin2z=0.0,
        delta_t=dt_data,
        f_lower=f_lower
    )
    t_gr = hp.sample_times.numpy()
    h_gr_raw = hp.numpy()

    peak_idx = np.argmax(np.abs(h_gr_raw))
    t_peak = t_gr[peak_idx]
    t_gr_shifted = t_gr - t_peak

    interp_gr = interp1d(t_gr_shifted, h_gr_raw, kind="cubic",
                         bounds_error=False, fill_value=0.0)
    out = interp_gr(t_data)
    return out


def pqg_bayesian_scan(y_fit, h_gr_fit, dt):
    """
    Vectorized PQG scan with Log-Likelihood & Marginal Likelihood (Evidence) integration.
    """
    N = len(h_gr_fit)
    H_gr = np.fft.rfft(h_gr_fit)
    freqs = np.fft.rfftfreq(N, d=dt)

    beta2_grid = np.linspace(beta2_min, beta2_max, beta2_N)
    b1_grid = np.linspace(b1_min, b1_max, N_b1)
    b0_grid = np.linspace(b0_min, b0_max, N_b0)

    amp = np.abs(H_gr)
    phi_gr = np.angle(H_gr)

    B1, B0 = np.meshgrid(b1_grid, b0_grid, indexing="ij")  # (N_b1, N_b0)
    B1_flat = B1.ravel()
    B0_flat = B0.ravel()

    rms_min_for_beta2 = {}
    params_for_beta2 = {}
    log_L_grid_b1_b0_per_beta2 = []

    freqs_2d = freqs[None, :]  # (1, N_freq)

    for i, beta2 in enumerate(beta2_grid):
        phi_beta2 = beta2 * freqs_2d**2  # (1, N_freq)
        total_phase = (phi_gr[None, :]
                       + phi_beta2
                       + B1_flat[:, None] * freqs_2d
                       + B0_flat[:, None])

        H_pqg_all = amp[None, :] * np.exp(1j * total_phase)
        h_pqg_all = np.fft.irfft(H_pqg_all, n=N, axis=1)

        diff = h_pqg_all - y_fit[None, :]
        rms_all = np.sqrt(np.mean(diff**2, axis=1))

        idx_min = np.argmin(rms_all)
        best_r = rms_all[idx_min]
        best_b1 = B1_flat[idx_min]
        best_b0 = B0_flat[idx_min]

        rms_min_for_beta2[beta2] = best_r
        params_for_beta2[beta2] = (best_b1, best_b0)

        # Reshape log_L to 2D grid (N_b1, N_b0)
        log_L_all = -0.5 * np.sum(diff**2, axis=1)
        log_L_grid_b1_b0_per_beta2.append(log_L_all.reshape((N_b1, N_b0)))

    return beta2_grid, b1_grid, b0_grid, rms_min_for_beta2, params_for_beta2, log_L_grid_b1_b0_per_beta2


# =========================================================
# 3) EVENT PROCESSING + BAYES FACTOR
# =========================================================

def process_event_with_approximant(args):
    name, t0, m1, m2, f_lower, fetch_win, approximant = args
    print(f"\n=== Processing event {name} with approximant {approximant} ===")

    try:
        data = {}
        for det in detectors:
            data[det] = load_and_condition(det, t0, fetch_win)

        dt = data["H1"]["dt"]

        h_gr_on_H1 = generate_gr_on_data_grid(approximant, dt, data["H1"]["t_data"], m1, m2, f_lower)
        h_gr_on_L1 = generate_gr_on_data_grid(approximant, dt, data["L1"]["t_data"], m1, m2, f_lower)

        results = {}

        for det, h_gr_on_data in zip(detectors, [h_gr_on_H1, h_gr_on_L1]):
            y_fit = data[det]["y_fit"]
            mask_fit = data[det]["mask_fit"]

            scale_GR = np.sqrt(np.mean(y_fit**2)) / (
                np.sqrt(np.mean(h_gr_on_data[mask_fit]**2)) + 1e-20
            )

            h_gr_fit_full = h_gr_on_data * scale_GR
            h_gr_fit = h_gr_fit_full[mask_fit]
            rms_GR = rms_error(y_fit, h_gr_fit)
            log_L_GR = compute_log_likelihood(y_fit, h_gr_fit)

            beta2_grid, b1_grid, b0_grid, rms_min_dict, params_dict, log_L_grids = pqg_bayesian_scan(y_fit, h_gr_fit, dt)

            results[det] = {
                "t_fit": data[det]["t_fit"],
                "y_fit": y_fit,
                "h_gr_fit": h_gr_fit,
                "rms_GR": rms_GR,
                "log_L_GR": log_L_GR,
                "beta2_grid": beta2_grid,
                "b1_grid": b1_grid,
                "b0_grid": b0_grid,
                "rms_min_dict": rms_min_dict,
                "params_dict": params_dict,
                "log_L_grids": log_L_grids,
            }

        # ---------------------------------------------------------
        # Combined Analysis & Bayes Factor Integration
        # ---------------------------------------------------------
        beta2_grid = results["H1"]["beta2_grid"]
        b1_grid = results["H1"]["b1_grid"]
        b0_grid = results["H1"]["b0_grid"]

        rms_H1_arr = np.array([results["H1"]["rms_min_dict"][b] for b in beta2_grid])
        rms_L1_arr = np.array([results["L1"]["rms_min_dict"][b] for b in beta2_grid])
        combined_rms = np.sqrt(0.5 * (rms_H1_arr**2 + rms_L1_arr**2))

        idx_best = np.argmin(combined_rms)
        beta2_best = beta2_grid[idx_best]

        H1_GR = results["H1"]["rms_GR"]
        L1_GR = results["L1"]["rms_GR"]
        comb_GR = np.sqrt(0.5 * (H1_GR**2 + L1_GR**2))
        comb_PQG = combined_rms[idx_best]

        # Total GR Log-Likelihood
        log_L_GR_total = results["H1"]["log_L_GR"] + results["L1"]["log_L_GR"]
        log_Z_GR = log_L_GR_total  # Delta prior for GR model

        # Numerical Evidence Integration for PQG
        prior_beta2 = 1.0 / (beta2_max - beta2_min)
        prior_b1 = 1.0 / (b1_max - b1_min)
        prior_b0 = 1.0 / (b0_max - b0_min)
        total_prior_density = prior_beta2 * (prior_b1**2) * (prior_b0**2)

        L_det_integrated = {"H1": np.zeros(beta2_N), "L1": np.zeros(beta2_N)}

        for det in detectors:
            for i in range(beta2_N):
                log_L_grid = results[det]["log_L_grids"][i]
                # Stabilized exponential for numerical integration
                L_grid_relative = np.exp(log_L_grid - results[det]["log_L_GR"])
                int_b0 = simpson(L_grid_relative, b0_grid, axis=1)
                int_b1_b0 = simpson(int_b0, b1_grid)
                L_det_integrated[det][i] = int_b1_b0

        # Combine detectors and integrate over beta2
        L_combined_beta2 = L_det_integrated["H1"] * L_det_integrated["L1"]
        Z_PQG_relative = simpson(L_combined_beta2, beta2_grid) * total_prior_density

        log_Z_PQG = log_L_GR_total + np.log(Z_PQG_relative + 1e-300)
        ln_B_PQG_GR = log_Z_PQG - log_Z_GR

        print(f"[EVENT {name}] Combined RMS: GR={comb_GR:.3e}, PQG={comb_PQG:.3e}")
        print(f"[EVENT {name}] Log evidence: GR={log_Z_GR:.2f}, PQG={log_Z_PQG:.2f} => ln B={ln_B_PQG_GR:.2f}")

        # PQG reconstruction for plotting
        for det in detectors:
            t_fit = results[det]["t_fit"]
            y_fit = results[det]["y_fit"]
            h_gr_fit = results[det]["h_gr_fit"]
            N_fit = len(h_gr_fit)
            H_gr = np.fft.rfft(h_gr_fit)
            freqs = np.fft.rfftfreq(N_fit, d=dt)
            b1_best, b0_best = results[det]["params_dict"][beta2_best]
            amp = np.abs(H_gr)
            phi_gr = np.angle(H_gr)
            total_phase = phi_gr + beta2_best * freqs**2 + b1_best * freqs + b0_best
            H_pqg = amp * np.exp(1j * total_phase)
            h_pqg_fit = np.fft.irfft(H_pqg, n=N_fit)
            rms_pqg = rms_error(y_fit, h_pqg_fit)
            results[det]["h_pqg_fit"] = h_pqg_fit
            results[det]["rms_PQG_best"] = rms_pqg

        # Time-domain plot: GR vs PQG-PN vs data
        plt.figure(figsize=(10, 6))
        for det, color_data in [("H1", "0.2"), ("L1", "0.6")]:
            t_fit = results[det]["t_fit"]
            y_fit = results[det]["y_fit"]
            h_gr_fit = results[det]["h_gr_fit"]
            h_pqg_fit = results[det]["h_pqg_fit"]

            plt.plot(t_fit, y_fit, ".", color=color_data, ms=2,
                     label=f"{det} data" if det == "H1" else None)
            plt.plot(t_fit, h_gr_fit, "orange", lw=1.3,
                     label="GR (fitted)" if det == "H1" else None)
            plt.plot(t_fit, h_pqg_fit, "blue", lw=1.3,
                     label=f"PQG-PN (beta2={beta2_best:.2e})" if det == "H1" else None)

        plt.xlabel("Time [s] relative to event center")
        plt.ylabel("Whitened strain [normalized]")
        plt.title(f"{name} ({approximant}): GR vs PQG-PN (ln B = {ln_B_PQG_GR:.2f})")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        fname_td = os.path.join(PLOT_DIR, f"{name}_{approximant}_GR_vs_PQGPN.png")
        plt.savefig(fname_td, dpi=200)
        plt.close()

        # RMS(beta2) plot
        plt.figure(figsize=(8, 5))
        plt.plot(beta2_grid, rms_H1_arr, "bo-", label="H1 RMS_min(beta2)")
        plt.plot(beta2_grid, rms_L1_arr, "go-", label="L1 RMS_min(beta2)")
        plt.plot(beta2_grid, combined_rms, "ko-", label="Combined RMS(beta2)")
        plt.axhline(H1_GR, color="red", ls="--", label="H1 RMS_GR")
        plt.axhline(L1_GR, color="orange", ls="--", label="L1 RMS_GR")
        plt.axhline(comb_GR, color="gold", ls="--", label="Combined RMS_GR")
        plt.axvline(beta2_best, color="black", ls=":", label=f"beta2_best={beta2_best:.2e}")
        plt.xlabel("beta2 (common PQG dispersion parameter)")
        plt.ylabel("RMS error")
        plt.title(f"{name} ({approximant}): RMS(beta2) [ln B={ln_B_PQG_GR:.2f}]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        fname_rms = os.path.join(PLOT_DIR, f"{name}_{approximant}_RMS_beta2.png")
        plt.savefig(fname_rms, dpi=200)
        plt.close()

        return {
            "name": name,
            "approximant": approximant,
            "beta2_best": beta2_best,
            "combined_RMS_GR": comb_GR,
            "combined_RMS_PQG": comb_PQG,
            "log_Z_GR": log_Z_GR,
            "log_Z_PQG": log_Z_PQG,
            "ln_B_PQG_GR": ln_B_PQG_GR,
        }

    except Exception as e:
        print(f"!!! Processing failed for {name}/{approximant}: {e}")
        return None


# =========================================================
# 4) FULL PIPELINE (PARALLELIZED)
# =========================================================

def run_rms_pipeline(EVENTS):
    print(f"[PIPELINE] Building task list from {len(EVENTS)} events")
    tasks = []
    for name, (t0, (m1, m2), f_lower, fetch_win) in EVENTS.items():
        for approx in APPROXIMANTS:
            tasks.append((name, t0, m1, m2, f_lower, fetch_win, approx))

    print(f"[PIPELINE] Total tasks: {len(tasks)}")
    print(f"[PIPELINE] CPU cores available: {cpu_count()}")

    with Pool(cpu_count()) as pool:
        results = pool.map(process_event_with_approximant, tasks)

    ok = [r for r in results if r is not None]
    print(f"[PIPELINE] Successful results: {len(ok)} / {len(tasks)}")
    return ok


# =========================================================
# 5) GLOBAL SIGNIFICANCE + CATALOG PLOTS
# =========================================================

def compute_global_significance(summaries):
    print(f"[STATS] Computing global significance from {len(summaries)} summaries")
    deltas = np.array([
        s["combined_RMS_GR"] - s["combined_RMS_PQG"]
        for s in summaries
    ])
    ln_B_vals = np.array([
        s["ln_B_PQG_GR"]
        for s in summaries
    ])

    mask = np.isfinite(deltas) & np.isfinite(ln_B_vals)
    deltas_clean = deltas[mask]
    ln_B_clean = ln_B_vals[mask]

    N = len(deltas_clean)
    print(f"[STATS] Finite data count: {N}")
    if N < 2:
        print("[STATS][WARN] Not enough samples for statistics")
        return {
            "N": N,
            "Z": np.nan,
            "p": np.nan,
            "mean_delta": np.nan,
            "std_delta": np.nan,
            "mean_ln_B": np.nan,
            "total_ln_B": np.nan,
            "deltas_clean": deltas_clean,
            "ln_B_clean": ln_B_clean,
        }

    mean_delta = np.mean(deltas_clean)
    std_delta = np.std(deltas_clean, ddof=1)
    stderr = std_delta / np.sqrt(N)

    Z = mean_delta / (stderr + 1e-30)
    p = 2 * (1 - 0.5 * (1 + erf(abs(Z) / sqrt(2))))

    mean_ln_B = np.mean(ln_B_clean)
    total_ln_B = np.sum(ln_B_clean)

    print(f"[STATS] mean_delta={mean_delta:.3e}, std_delta={std_delta:.3e}, Z={Z:.3f}, p={p:.3e}")
    print(f"[STATS] mean ln B={mean_ln_B:.2f}, total cumulative ln B={total_ln_B:.2f}")

    return {
        "N": N,
        "mean_delta": mean_delta,
        "std_delta": std_delta,
        "Z": Z,
        "p": p,
        "mean_ln_B": mean_ln_B,
        "total_ln_B": total_ln_B,
        "deltas_clean": deltas_clean,
        "ln_B_clean": ln_B_clean,
    }


def write_summary(summaries, stats, filename="summary.txt"):
    print(f"[IO] Writing summary file: {filename}")
    with open(filename, "w") as f:
        f.write("PQG-PN vs GR RMS & Bayes Factor summary\n")
        f.write(f"Generated: {datetime.now()}\n\n")

        for s in summaries:
            f.write(f"Event: {s['name']}\n")
            f.write(f"  Approximant: {s['approximant']}\n")
            f.write(f"  beta2_best: {s['beta2_best']:.6e}\n")
            f.write(f"  Combined RMS GR:  {s['combined_RMS_GR']:.6e}\n")
            f.write(f"  Combined RMS PQG: {s['combined_RMS_PQG']:.6e}\n")
            f.write(f"  Delta RMS: {s['combined_RMS_GR'] - s['combined_RMS_PQG']:.6e}\n")
            f.write(f"  Log Evidence GR:  {s['log_Z_GR']:.4f}\n")
            f.write(f"  Log Evidence PQG: {s['log_Z_PQG']:.4f}\n")
            f.write(f"  ln B (PQG/GR):    {s['ln_B_PQG_GR']:.4f}\n\n")

        f.write("\nGLOBAL SIGNIFICANCE:\n")
        f.write(f"  N pairs: {stats['N']}\n")
        f.write(f"  mean Δ RMS: {stats['mean_delta']:.6e}\n")
        f.write(f"  std Δ RMS:  {stats['std_delta']:.6e}\n")
        f.write(f"  Z-score: {stats['Z']:.3f}\n")
        f.write(f"  p-value: {stats['p']:.3e}\n")
        f.write(f"  mean ln B (PQG/GR): {stats['mean_ln_B']:.4f}\n")
        f.write(f"  total cumulative ln B: {stats['total_ln_B']:.4f}\n")


def write_summary_clean(stats, filename="summary_clean.txt"):
    print(f"[IO] Writing cleaned summary file: {filename}")
    with open(filename, "w") as f:
        f.write("GLOBAL SIGNIFICANCE (cleaned, NaN-free with Bayes Factor)\n\n")
        f.write(f"N pairs: {stats['N']}\n")
        f.write(f"mean Δ RMS: {stats['mean_delta']:.6e}\n")
        f.write(f"std Δ RMS:  {stats['std_delta']:.6e}\n")
        f.write(f"Z-score: {stats['Z']:.3f}\n")
        f.write(f"p-value: {stats['p']:.3e}\n")
        f.write(f"mean ln B (PQG/GR): {stats['mean_ln_B']:.4f}\n")
        f.write(f"total cumulative ln B: {stats['total_ln_B']:.4f}\n")


def make_catalog_plots(summaries, stats):
    print("[PLOT] Creating catalog-level plots")
    # --- beta2 histogram per approximant
    by_approx = {}
    for s in summaries:
        approx = s["approximant"]
        by_approx.setdefault(approx, []).append(s["beta2_best"])

    for approx, beta_list in by_approx.items():
        beta_arr = np.array(beta_list)
        med = np.median(beta_arr)
        lo = np.percentile(beta_arr, 16)
        hi = np.percentile(beta_arr, 84)

        plt.figure(figsize=(7, 4))
        plt.hist(beta_arr, bins=40, density=True, color="C0", alpha=0.7)
        plt.axvline(med, color="k", ls="-", label="median")
        plt.axvline(lo, color="k", ls="--", label="16/84% CI")
        plt.axvline(hi, color="k", ls="--")
        plt.xlabel(r"$\beta_2$ [s$^2$]")
        plt.ylabel("posterior density (empirical)")
        plt.title(rf"Catalog-level $\beta_2$ distribution ({approx})")
        plt.legend()
        plt.tight_layout()
        fname = os.path.join(PLOT_DIR, f"beta2_posterior_{approx}.png")
        plt.savefig(fname, dpi=200)
        plt.close()

    # --- Δ distribution + normal fit
    deltas = stats["deltas_clean"]
    if len(deltas) > 0:
        mu = np.mean(deltas)
        sigma = np.std(deltas, ddof=1)

        plt.figure(figsize=(7, 4))
        plt.hist(deltas, bins=40, density=True, alpha=0.7, color="C0", label="Δ RMS histogram")
        x = np.linspace(mu - 4*sigma, mu + 4*sigma, 400)
        gauss = (1.0 / (sigma * np.sqrt(2*np.pi))) * np.exp(-0.5 * ((x - mu)/sigma)**2)
        plt.plot(x, gauss, "r-", lw=2, label="Normal fit")
        plt.xlabel(r"$\Delta = \mathrm{RMS}_{GR} - \mathrm{RMS}_{PQG}$")
        plt.ylabel("density")
        plt.title("Distribution of Δ RMS across catalog")
        plt.legend()
        plt.tight_layout()
        fname = os.path.join(PLOT_DIR, "delta_distribution.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        # --- Bayes Factor (ln B) distribution plot
        ln_B_vals = stats["ln_B_clean"]
        plt.figure(figsize=(7, 4))
        plt.hist(ln_B_vals, bins=35, color="teal", alpha=0.7, edgecolor="black")
        plt.axvline(0, color="red", linestyle="--", label="ln B = 0 (GR baseline)")
        plt.axvline(np.mean(ln_B_vals), color="black", linestyle="-", label=f"Mean ln B = {np.mean(ln_B_vals):.2f}")
        plt.xlabel(r"$\ln \mathcal{B}_{\mathrm{PQG/GR}}$ (Log Bayes Factor)")
        plt.ylabel("Event count")
        plt.title(r"Catalog Distribution of Log Bayes Factors ($\ln \mathcal{B}_{\mathrm{PQG/GR}}$)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        fname_bf = os.path.join(PLOT_DIR, "bayes_factor_distribution.png")
        plt.savefig(fname_bf, dpi=200)
        plt.close()

        # --- cumulative Δ(N)
        N = len(deltas)
        cum_mean = np.array([np.mean(deltas[:i+1]) for i in range(N)])

        plt.figure(figsize=(7, 4))
        plt.plot(np.arange(1, N+1), cum_mean, "k-")
        plt.xlabel("Number of event–approximant pairs N")
        plt.ylabel(r"Cumulative mean $\Delta(N)$")
        plt.title("Cumulative improvement Δ(N)")
        plt.grid(True)
        plt.tight_layout()
        fname = os.path.join(PLOT_DIR, "cumulative_delta.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        # --- Z-score growth Z(N)
        Z_vals = []
        for i in range(N):
            if i == 0:
                Z_vals.append(0.0)
            else:
                sub = deltas[:i+1]
                mu_i = np.mean(sub)
                sigma_i = np.std(sub, ddof=1)
                stderr_i = sigma_i / np.sqrt(len(sub))
                Z_i = mu_i / (stderr_i + 1e-30)
                Z_vals.append(Z_i)
        Z_vals = np.array(Z_vals)

        plt.figure(figsize=(7, 4))
        plt.plot(np.arange(1, N+1), Z_vals, "k-")
        plt.axhline(stats["Z"], color="r", ls="--", label=f"Final Z = {stats['Z']:.3f}")
        plt.xlabel("Number of event–approximant pairs N")
        plt.ylabel("Z(N)")
        plt.title("Growth of Z-score across catalog")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        fname = os.path.join(PLOT_DIR, "zscore_growth.png")
        plt.savefig(fname, dpi=200)
        plt.close()


# =========================================================
# 6) MAIN
# =========================================================

def main():
    print("=== FULL_CATALOG_SCAN_PLOT WITH BAYES FACTOR START ===")

    EVENTS = load_events_from_csv("events.csv")
    print(f"[MAIN] Number of loaded events: {len(EVENTS)}")

    summaries = run_rms_pipeline(EVENTS)
    print(f"[MAIN] Summaries collected: {len(summaries)}")

    stats = compute_global_significance(summaries)

    write_summary(summaries, stats, filename="summary.txt")
    write_summary_clean(stats, filename="summary_clean.txt")
    make_catalog_plots(summaries, stats)

    print("=== FULL_CATALOG_SCAN_PLOT FINISHED ===")
    print("summary.txt, summary_clean.txt and PNG plots are available in the 'plots_catalog' directory.")


if __name__ == "__main__":
    main()