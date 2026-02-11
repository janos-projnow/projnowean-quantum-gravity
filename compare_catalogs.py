# -*- coding: utf-8 -*-
"""
Comparative statistical analysis script:
  - summary.txt       (real GWOSC catalog)
  - summary_null.txt  (null test: GR + noise)
Tasks:
  * Δ-histogram comparison
  * β2 distribution comparison
  * computation of global Z-scores
  * cumulative Δ(N) curves
"""

import numpy as np
import matplotlib.pyplot as plt
import re
from math import erf, sqrt


REAL_SUMMARY = "summary.txt"
NULL_SUMMARY = "summary_null.txt"


# ---------------------------------------------------------
# 1) Summary parser
# ---------------------------------------------------------

def parse_summary(path):
    beta2_list = []
    delta_list = []

    beta2_re = re.compile(r"beta2_best:\s*([+-]?\d+\.\d+e[+-]\d+)")
    delta_re = re.compile(r"Delta:\s*([+-]?\d+\.\d+e[+-]\d+|nan)")

    with open(path, "r") as f:
        lines = f.readlines()

    current_beta2 = None
    for line in lines:
        m_b = beta2_re.search(line)
        if m_b:
            try:
                current_beta2 = float(m_b.group(1))
            except:
                current_beta2 = None
            continue

        m_d = delta_re.search(line)
        if m_d and current_beta2 is not None:
            try:
                d = float(m_d.group(1))
            except:
                d = np.nan

            if np.isfinite(d) and np.isfinite(current_beta2):
                beta2_list.append(current_beta2)
                delta_list.append(d)

            current_beta2 = None

    return np.array(beta2_list), np.array(delta_list)


# ---------------------------------------------------------
# 2) Global significance (same formula as in the pipeline)
# ---------------------------------------------------------

def compute_global_significance(deltas):
    deltas = np.asarray(deltas)
    deltas = deltas[np.isfinite(deltas)]

    N = len(deltas)
    if N < 2:
        return {"N": N, "mean": np.nan, "std": np.nan, "Z": np.nan, "p": np.nan}

    mean = np.mean(deltas)
    std = np.std(deltas, ddof=1)
    stderr = std / np.sqrt(N)

    Z = mean / (stderr + 1e-30)
    p = 2 * (1 - 0.5 * (1 + erf(abs(Z) / sqrt(2))))

    return {"N": N, "mean": mean, "std": std, "Z": Z, "p": p}


# ---------------------------------------------------------
# 3) Plot cumulative Δ(N)
# ---------------------------------------------------------

def plot_cumulative(delta_real, delta_null):
    plt.figure(figsize=(7,5))

    # Real catalog
    cum_real = np.cumsum(delta_real) / np.arange(1, len(delta_real)+1)
    plt.plot(cum_real, label="Real catalog", color="C1")

    # Null test
    cum_null = np.cumsum(delta_null) / np.arange(1, len(delta_null)+1)
    plt.plot(cum_null, label="Null test", color="C0")

    plt.xlabel("N (event–template pairs)")
    plt.ylabel("Cumulative mean Δ(N)")
    plt.title("Cumulative Δ(N): real vs null test")
    plt.legend()
    plt.tight_layout()
    plt.savefig("cumulative_delta_compare.png", dpi=200)


# ---------------------------------------------------------
# 4) Main
# ---------------------------------------------------------

def main():
    print("=== Comparative statistical analysis ===")

    beta2_real, delta_real = parse_summary(REAL_SUMMARY)
    beta2_null, delta_null = parse_summary(NULL_SUMMARY)

    # --- Global statistics ---
    stats_real = compute_global_significance(delta_real)
    stats_null = compute_global_significance(delta_null)

    print("\nREAL CATALOG:")
    for k, v in stats_real.items():
        print(f"  {k}: {v}")

    print("\nNULL TEST:")
    for k, v in stats_null.items():
        print(f"  {k}: {v}")

    # --- Δ histogram ---
    plt.figure(figsize=(7,5))
    plt.hist(delta_null, bins=40, density=True, alpha=0.5, label="Null test", color="C0")
    plt.hist(delta_real, bins=40, density=True, alpha=0.5, label="Real catalog", color="C1")
    plt.xlabel(r"$\Delta = \mathrm{RMS}_{GR} - \mathrm{RMS}_{PQG}$")
    plt.ylabel("Density")
    plt.title("Δ distribution: real vs null test")
    plt.legend()
    plt.tight_layout()
    plt.savefig("delta_hist_compare.png", dpi=200)

    # --- β2 histogram ---
    plt.figure(figsize=(7,5))
    plt.hist(beta2_null, bins=40, density=True, alpha=0.5,
             label=r"$\beta_2$ null test", color="C0")
    plt.hist(beta2_real, bins=40, density=True, alpha=0.5,
             label=r"$\beta_2$ real catalog", color="C1")
    plt.xlabel(r"$\beta_2\ [s^2]$")
    plt.ylabel("Density")
    plt.title(r"$\beta_2$ distribution: real vs null test")
    plt.legend()
    plt.tight_layout()
    plt.savefig("beta2_hist_compare.png", dpi=200)

    # --- Cumulative Δ(N) ---
    plot_cumulative(delta_real, delta_null)

    print("\nFigures saved:")
    print("  delta_hist_compare.png")
    print("  beta2_hist_compare.png")
    print("  cumulative_delta_compare.png")


if __name__ == "__main__":
    main()
