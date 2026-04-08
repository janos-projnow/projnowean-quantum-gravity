#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import math
import os
from datetime import datetime, UTC

# ---------------------------------------------
# Output directory
# ---------------------------------------------
OUTDIR = "PQG_null_test_results"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------
# PQG model (same as photon pipeline)
# ---------------------------------------------

def pqg_time_delay(beta2, energies, ref_energy):
    return beta2 * (energies - ref_energy)**2

def model_emission_times(beta2, times_rel, energies, ref_energy):
    return times_rel - pqg_time_delay(beta2, energies, ref_energy)

def log_likelihood_beta2(beta2, times_rel, energies, ref_energy):
    if len(times_rel) == 0:
        return -np.inf
    t_emit = model_emission_times(beta2, times_rel, energies, ref_energy)
    t_centered = t_emit - np.mean(t_emit)
    var = np.var(t_centered)
    sigma_model = 0.1
    return -0.5 * var / (sigma_model**2)

def log_prior_beta2(beta2, mean, sigma):
    return -0.5 * ((beta2 - mean)**2) / (sigma**2)

# ---------------------------------------------
# Posterior grid (same as photon8)
# ---------------------------------------------

def compute_posterior_grid(times_rel, energies,
                           beta2_prior_mean=0.0,
                           beta2_prior_sigma=1e-4,
                           beta2_min=-5e-4,
                           beta2_max=5e-4,
                           n_grid=1001):

    ref_energy = np.median(energies)
    beta2_grid = np.linspace(beta2_min, beta2_max, n_grid)
    log_post = []

    for b in tqdm(beta2_grid, desc="β2 grid"):
        lp = log_prior_beta2(b, beta2_prior_mean, beta2_prior_sigma)
        ll = log_likelihood_beta2(b, times_rel, energies, ref_energy)
        log_post.append(lp + ll)

    log_post = np.array(log_post)
    log_post -= np.max(log_post)
    post = np.exp(log_post)
    post /= np.trapezoid(post, beta2_grid)

    mean = np.trapezoid(beta2_grid * post, beta2_grid)
    var = np.trapezoid((beta2_grid - mean)**2 * post, beta2_grid)
    sigma = math.sqrt(var)

    return beta2_grid, post, mean, sigma

# ---------------------------------------------
# Synthetic noise generator
# ---------------------------------------------

def generate_synthetic_noise(N=300000,
                             time_range=(0, 50),
                             energy_range=(0, 120),
                             seed=None):

    rng = np.random.default_rng(seed)

    # Gaussian or uniform time noise
    times = rng.normal(loc=25, scale=10, size=N)
    times = np.clip(times, time_range[0], time_range[1])

    # Uniform energy noise
    energies = rng.uniform(energy_range[0], energy_range[1], size=N)

    return times, energies

# ---------------------------------------------
# Main null test
# ---------------------------------------------

def run_null_test(N=300000, seed=1234):

    print(f"\n=== Running PQG NULL TEST (N={N}, seed={seed}) ===")

    times_rel, energies = generate_synthetic_noise(N=N, seed=seed)

    beta2_grid, post, mean, sigma = compute_posterior_grid(
        times_rel, energies,
        beta2_prior_mean=0.0,
        beta2_prior_sigma=1e-4,
        beta2_min=-5e-4,
        beta2_max=5e-4,
        n_grid=1001,
    )

    Z = mean / sigma

    print(f"  β2_mean  = {mean:.3e}")
    print(f"  β2_sigma = {sigma:.3e}")
    print(f"  Z-score  = {Z:.3f}")

    # Save posterior plot
    plt.figure(figsize=(6,4))
    plt.plot(beta2_grid, post)
    plt.axvline(mean, color="orange", linestyle="--")
    plt.title(f"NULL TEST posterior (N={N})")
    plt.xlabel("β2")
    plt.ylabel("p(β2|noise)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f"null_posterior_N{N}.png"), dpi=200)
    plt.close()

    return mean, sigma, Z

# ---------------------------------------------
# Run multiple seeds
# ---------------------------------------------

def main():
    results = []
    for seed in [1,2,3,4,5,6,7,8,9,10]:
        mean, sigma, Z = run_null_test(N=300000, seed=seed)
        results.append(Z)

    print("\n=== NULL TEST SUMMARY ===")
    print("Z-scores:", results)
    print(f"Mean Z = {np.mean(results):.3f}")
    print(f"Std  Z = {np.std(results):.3f}")

if __name__ == "__main__":
    main()
