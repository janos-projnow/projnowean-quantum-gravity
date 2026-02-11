# -*- coding: utf-8 -*-
"""
PQG-PN MASTER SCRIPT – CSV-based catalog analysis, cache + vectorized PQG scan

Functions:
  - Loads the GW event list from events.csv
  - For each event and approximant runs:
        GR RMS  vs  PQG-PN RMS(beta2)
  - TimeSeries.fetch_open_data cached (by detector, t0, fetch_win)
  - pqg_scan vectorized (b0-b1 grid at once, loop only over beta2)
  - Does not generate plots
  - In summary.txt:
        * per-event results
        * global RMS improvement significance (Z-score, p-value)
"""

import numpy as np
import csv
from datetime import datetime
from gwpy.timeseries import TimeSeries
from pycbc.waveform import get_td_waveform
from scipy.interpolate import interp1d
from math import erf, sqrt
from multiprocessing import Pool, cpu_count

# =========================================================
# 0) PARAMETERS
# =========================================================

APPROXIMANTS = [
    "IMRPhenomD",
    "SEOBNRv4",
]

detectors = ["H1", "L1"]

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

# fetch_open_data cache (det, t_start, t_end) -> TimeSeries
_FETCH_CACHE = {}


# =========================================================
# 1) CSV LOADING
# =========================================================

def load_events_from_csv(csv_path="events.csv", max_events=None):
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
            except Exception:
                continue

            f_lower = 30.0
            fetch_win = 4.0

            EVENTS[name] = (gps, (m1, m2), f_lower, fetch_win)
            count += 1

            if max_events is not None and count >= max_events:
                break

    return EVENTS


# =========================================================
# 2) HELPER FUNCTIONS
# =========================================================

def rms_error(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))


def fetch_open_data_cached(det, t_start, t_end):
    key = (det, float(t_start), float(t_end))
    if key in _FETCH_CACHE:
        return _FETCH_CACHE[key]
    h = TimeSeries.fetch_open_data(det, t_start, t_end, cache=True)
    _FETCH_CACHE[key] = h
    return h


def load_and_condition(det, t0, fetch_win):
    t_start = t0 - fetch_win/2.0
    t_end   = t0 + fetch_win/2.0
    print(f"[{det}] Fetching {t_start} .. {t_end}")

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
    return interp_gr(t_data)


def pqg_scan(y_fit, h_gr_fit, dt):
    """
    Vectorized PQG scan:
      - Single FFT
      - b0-b1 grid vectorized
      - Loop only over beta2
    """
    N = len(h_gr_fit)
    H_gr = np.fft.rfft(h_gr_fit)
    freqs = np.fft.rfftfreq(N, d=dt)

    beta2_grid = np.linspace(beta2_min, beta2_max, beta2_N)
    b1_grid = np.linspace(b1_min, b1_max, N_b1)
    b0_grid = np.linspace(b0_min, b0_max, N_b0)

    amp = np.abs(H_gr)
    phi_gr = np.angle(H_gr)

    # b0-b1 grid
    B1, B0 = np.meshgrid(b1_grid, b0_grid, indexing="ij")  # shape (N_b1, N_b0)
    B1_flat = B1.ravel()  # (N_b1*N_b0,)
    B0_flat = B0.ravel()

    rms_min_for_beta2 = {}
    params_for_beta2 = {}

    # freqs shape: (N_freq,)
    # B1_flat, B0_flat shape: (N_phase,)
    # -> broadcast: (N_phase, N_freq)
    freqs_2d = freqs[None, :]  # (1, N_freq)

    for beta2 in beta2_grid:
        phi_beta2 = beta2 * freqs_2d**2  # (1, N_freq)
        # total phase for all (b1,b0) combinations:
        # phi_gr: (N_freq,) -> (1, N_freq)
        total_phase = (phi_gr[None, :]
                       + phi_beta2
                       + B1_flat[:, None] * freqs_2d
                       + B0_flat[:, None])

        H_pqg_all = amp[None, :] * np.exp(1j * total_phase)  # (N_phase, N_freq)
        h_pqg_all = np.fft.irfft(H_pqg_all, n=N, axis=1)     # (N_phase, N_time)

        # RMS for each (b1,b0) combination
        diff = h_pqg_all - y_fit[None, :]                    # (N_phase, N_time)
        rms_all = np.sqrt(np.mean(diff**2, axis=1))          # (N_phase,)

        idx_min = np.argmin(rms_all)
        best_r = rms_all[idx_min]
        best_b1 = B1_flat[idx_min]
        best_b0 = B0_flat[idx_min]

        rms_min_for_beta2[beta2] = best_r
        params_for_beta2[beta2] = (best_b1, best_b0)

    return beta2_grid, rms_min_for_beta2, params_for_beta2


# =========================================================
# 3) EVENT PROCESSING
# =========================================================

def process_event_with_approximant(args):
    name, t0, m1, m2, f_lower, fetch_win, approximant = args
    print(f"\n=== {name} / {approximant} ===")

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

            beta2_grid, rms_min_dict, params_dict = pqg_scan(y_fit, h_gr_fit, dt)

            results[det] = {
                "rms_GR": rms_GR,
                "beta2_grid": beta2_grid,
                "rms_min_dict": rms_min_dict,
                "params_dict": params_dict,
            }

        beta2_grid = results["H1"]["beta2_grid"]
        rms_H1_arr = np.array([results["H1"]["rms_min_dict"][b] for b in beta2_grid])
        rms_L1_arr = np.array([results["L1"]["rms_min_dict"][b] for b in beta2_grid])
        combined_rms = np.sqrt(0.5 * (rms_H1_arr**2 + rms_L1_arr**2))

        idx_best = np.argmin(combined_rms)
        beta2_best = beta2_grid[idx_best]

        H1_GR = results["H1"]["rms_GR"]
        L1_GR = results["L1"]["rms_GR"]
        comb_GR = np.sqrt(0.5 * (H1_GR**2 + L1_GR**2))
        comb_PQG = combined_rms[idx_best]

        return {
            "name": name,
            "approximant": approximant,
            "beta2_best": beta2_best,
            "combined_RMS_GR": comb_GR,
            "combined_RMS_PQG": comb_PQG,
        }

    except Exception as e:
        print(f"!!! Failed {name}/{approximant}: {e}")
        return None


# =========================================================
# 4) FULL PIPELINE (PARALLELIZED)
# =========================================================

def run_rms_pipeline(EVENTS):
    tasks = []
    for name, (t0, (m1, m2), f_lower, fetch_win) in EVENTS.items():
        for approx in APPROXIMANTS:
            tasks.append((name, t0, m1, m2, f_lower, fetch_win, approx))

    print(f"[INFO] Total tasks: {len(tasks)}")
    print(f"[INFO] CPU cores: {cpu_count()}")

    with Pool(cpu_count()) as pool:
        results = pool.map(process_event_with_approximant, tasks)

    return [r for r in results if r is not None]


# =========================================================
# 5) GLOBAL SIGNIFICANCE
# =========================================================

def compute_global_significance(summaries):
    deltas = np.array([
        s["combined_RMS_GR"] - s["combined_RMS_PQG"]
        for s in summaries
    ])

    N = len(deltas)
    if N < 2:
        return {"N": N, "Z": np.nan, "p": np.nan,
                "mean_delta": np.nan, "std_delta": np.nan}

    mean_delta = np.mean(deltas)
    std_delta = np.std(deltas, ddof=1)
    stderr = std_delta / np.sqrt(N)

    Z = mean_delta / (stderr + 1e-30)
    p = 2 * (1 - 0.5 * (1 + erf(abs(Z) / sqrt(2))))

    return {
        "N": N,
        "mean_delta": mean_delta,
        "std_delta": std_delta,
        "Z": Z,
        "p": p,
    }


# =========================================================
# 6) SUMMARY WRITING
# =========================================================

def write_summary(summaries, stats, filename="summary.txt"):
    with open(filename, "w") as f:
        f.write("PQG-PN vs GR RMS summary\n")
        f.write(f"Generated: {datetime.now()}\n\n")

        for s in summaries:
            f.write(f"Event: {s['name']}\n")
            f.write(f"  Approximant: {s['approximant']}\n")
            f.write(f"  beta2_best: {s['beta2_best']:.6e}\n")
            f.write(f"  Combined RMS GR:  {s['combined_RMS_GR']:.6e}\n")
            f.write(f"  Combined RMS PQG: {s['combined_RMS_PQG']:.6e}\n")
            f.write(f"  Delta: {s['combined_RMS_GR'] - s['combined_RMS_PQG']:.6e}\n\n")

        f.write("\nGLOBAL SIGNIFICANCE:\n")
        f.write(f"  N pairs: {stats['N']}\n")
        f.write(f"  mean Δ: {stats['mean_delta']:.6e}\n")
        f.write(f"  std Δ:  {stats['std_delta']:.6e}\n")
        f.write(f"  Z-score: {stats['Z']:.3f}\n")
        f.write(f"  p-value: {stats['p']:.3e}\n")


# =========================================================
# 7) MAIN
# =========================================================

def main():
    print("=== PQG MASTER PIPELINE START ===")

    EVENTS = load_events_from_csv("events.csv")
    print(f"[INFO] Loaded events: {len(EVENTS)}")

    summaries = run_rms_pipeline(EVENTS)
    stats = compute_global_significance(summaries)
    write_summary(summaries, stats)

    print("=== PQG MASTER PIPELINE DONE ===")
    print("summary.txt generated.")


if __name__ == "__main__":
    main()

