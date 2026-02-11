# -*- coding: utf-8 -*-
# null_test.py
# Null simulation for the PQG-PN pipeline: pure GR signal + Gaussian noise
# 300 synthetic events, IMRPhenomD and SEOBNRv4, with GR injection

import numpy as np
from datetime import datetime
from pycbc.waveform import get_td_waveform
from scipy.interpolate import interp1d
from math import erf, sqrt

# ====== PARAMETERS ======

N_EVENTS = 300

APPROXIMANTS = [
    "IMRPhenomD",
    "SEOBNRv4",
]

whiten_seg = 4.0      # kept for consistency; no actual whitening is performed
f_low = 30.0
f_high = 300.0        # not explicitly used, retained for structural consistency

fit_tmin = -0.10
fit_tmax =  0.02

beta2_min = -1.0e-4
beta2_max =  1.0e-4
beta2_N   = 41

N_b1 = 1   # b1 fixed to 0 (only β2 is free)
N_b0 = 1   # b0 fixed to 0

b1_min = 0.0
b1_max = 0.0

b0_min = 0.0
b0_max = 0.0

# time grid parameters
dt = 1.0 / 4096.0
t_window = 0.25   # total window length [s]
t0 = 0.0          # signal peak around t = 0

noise_sigma = 1.0  # Gaussian noise standard deviation (relative scale; sets SNR)

rng = np.random.default_rng(12345)


# ====== HELPER FUNCTIONS ======

def rms_error(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))


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

    interp_gr = interp1d(
        t_gr_shifted, h_gr_raw, kind="cubic",
        bounds_error=False, fill_value=0.0
    )
    return interp_gr(t_data)


def pqg_scan(y_fit, h_gr_fit, dt):
    """
    PQG scan for a pure-GR null test:
    - only beta2 is free (b0 = b1 = 0)
    """
    N = len(h_gr_fit)
    H_gr = np.fft.rfft(h_gr_fit)
    freqs = np.fft.rfftfreq(N, d=dt)

    beta2_grid = np.linspace(beta2_min, beta2_max, beta2_N)
    b1_grid = np.linspace(b1_min, b1_max, N_b1)
    b0_grid = np.linspace(b0_min, b0_max, N_b0)

    amp = np.abs(H_gr)
    phi_gr = np.angle(H_gr)

    B1, B0 = np.meshgrid(b1_grid, b0_grid, indexing="ij")
    B1_flat = B1.ravel()
    B0_flat = B0.ravel()

    rms_min_for_beta2 = {}
    params_for_beta2 = {}

    freqs_2d = freqs[None, :]

    for beta2 in beta2_grid:
        phi_beta2 = beta2 * freqs_2d**2
        total_phase = (
            phi_gr[None, :]
            + phi_beta2
            + B1_flat[:, None] * freqs_2d
            + B0_flat[:, None]
        )

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

    return beta2_grid, rms_min_for_beta2, params_for_beta2


def compute_global_significance(summaries):
    deltas = np.array([
        s["combined_RMS_GR"] - s["combined_RMS_PQG"]
        for s in summaries
    ])

    N = len(deltas)
    if N < 2:
        return {
            "N": N,
            "Z": np.nan,
            "p": np.nan,
            "mean_delta": np.nan,
            "std_delta": np.nan,
        }

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


# ====== MAIN NULL-TEST PIPELINE ======

def main():
    print("=== PQG NULL TEST START ===")

    summaries = []

    # time grid
    N_samples = int(t_window / dt)
    t_data = np.linspace(-t_window/2, t_window/2, N_samples, endpoint=False)

    mask_fit = (t_data > fit_tmin) & (t_data < fit_tmax)
    t_fit = t_data[mask_fit]

    for i in range(N_EVENTS):
        # mixture of symmetric / asymmetric BH binaries
        if i % 2 == 0:
            # symmetric: q ~ 1
            m1 = rng.uniform(20, 60)
            q = rng.uniform(0.8, 1.0)
        else:
            # asymmetric: q << 1
            m1 = rng.uniform(20, 60)
            q = rng.uniform(0.2, 0.7)
        m2 = m1 * q

        for approximant in APPROXIMANTS:
            name = f"SIM_{i:03d}_{approximant}"
            print(f"\n=== {name} ===")

            try:
                # GR signal generation (injection)
                h_gr = generate_gr_on_data_grid(
                    approximant, dt, t_data, m1, m2, f_low
                )

                # add Gaussian noise
                noise = rng.normal(0.0, noise_sigma, size=len(t_data))
                y_data = h_gr + noise

                # fit window
                y_fit = y_data[mask_fit]
                h_gr_fit = h_gr[mask_fit]

                # amplitude rescaling (same as in the original code)
                scale_GR = np.sqrt(np.mean(y_fit**2)) / (
                    np.sqrt(np.mean(h_gr_fit**2)) + 1e-20
                )
                h_gr_fit_full = h_gr * scale_GR
                h_gr_fit = h_gr_fit_full[mask_fit]

                rms_GR = rms_error(y_fit, h_gr_fit)

                # PQG scan (only beta2 free)
                beta2_grid, rms_min_dict, params_dict = pqg_scan(
                    y_fit, h_gr_fit, dt
                )

                rms_arr = np.array([rms_min_dict[b] for b in beta2_grid])
                idx_best = np.argmin(rms_arr)
                beta2_best = beta2_grid[idx_best]
                rms_PQG = rms_arr[idx_best]

                summaries.append({
                    "name": name,
                    "approximant": approximant,
                    "beta2_best": beta2_best,
                    "combined_RMS_GR": rms_GR,
                    "combined_RMS_PQG": rms_PQG,
                })

            except Exception as e:
                print(f"!!! Failed {name}: {e}")
                continue

    stats = compute_global_significance(summaries)

    # write summary_null.txt
    with open("summary_null.txt", "w") as f:
        f.write("PQG-PN vs GR RMS NULL TEST summary (synthetic GR + noise)\n")
        f.write(f"Generated: {datetime.now()}\n\n")

        for s in summaries:
            f.write(f"Event: {s['name']}\n")
            f.write(f"  Approximant: {s['approximant']}\n")
            f.write(f"  beta2_best: {s['beta2_best']:.6e}\n")
            f.write(f"  Combined RMS GR:  {s['combined_RMS_GR']:.6e}\n")
            f.write(f"  Combined RMS PQG: {s['combined_RMS_PQG']:.6e}\n")
            f.write(
                f"  Delta: {s['combined_RMS_GR'] - s['combined_RMS_PQG']:.6e}\n\n"
            )

        f.write("\nGLOBAL SIGNIFICANCE (NULL TEST):\n")
        f.write(f"  N pairs: {stats['N']}\n")
        f.write(f"  mean Δ: {stats['mean_delta']:.6e}\n")
        f.write(f"  std Δ:  {stats['std_delta']:.6e}\n")
        f.write(f"  Z-score: {stats['Z']:.3f}\n")
        f.write(f"  p-value: {stats['p']:.3e}\n")

    print("=== PQG NULL TEST DONE ===")
    print("summary_null.txt generated.")


if __name__ == "__main__":
    main()
