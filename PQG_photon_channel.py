#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PQG photon-channel pipeline for multi-event GRB analysis.

This script implements a complete, self-contained analysis chain that:
  - downloads Fermi/GBM TTE data for a catalog of GRBs,
  - infers a trigger time and analysis window (a simple "mini MET-seeker"),
  - loads photon arrival times and energies within that window,
  - fits a quadratic Planck-scale quantum gravity (PQG) time-delay parameter β₂
    using a simple likelihood model on a 1D grid with adaptive range,
  - produces posterior plots and time-of-arrival vs energy diagnostic plots,
  - writes per-event text summaries,
  - and combines all events into a single inverse-variance–weighted β₂ constraint,
    reporting a global Z-score (β₂ / σ_β₂).

The code is written to be readable for physicists who may not be familiar
with the internal development history of the pipeline.
"""

import os
import csv
import math
import urllib.request
from datetime import datetime, UTC

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm

# -------------------------------------------------------------------------
# Global output directory for all photon-channel results
# -------------------------------------------------------------------------
BASE_OUTPUT_DIR = "PQG_Photon_channel"
os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

# -------------------------------------------------------------------------
# 0) Event catalog
# -------------------------------------------------------------------------
# Each entry defines:
#   - a human-readable GRB name,
#   - the corresponding Fermi/GBM trigger ID,
#   - the year (used in the FSSC directory structure),
#   - the list of GBM detectors for which we attempt to download TTE data,
#   - optional manual trigger MET and time window (if known a priori).
#
# If manual_trigger_met and manual_time_window are None, the pipeline will
# estimate them from the TTE data using a simple "mini MET-seeker".
EVENT_CATALOG = [
    {
        "name": "GRB170817A",
        "gbm_trigger": "bn170817529",
        "year": "2017",
        "tte_detectors": ["n0", "n1"],
        "manual_trigger_met": 524666403.0,
        "manual_time_window": (-5.0, 15.0),
    },
    {
        "name": "GRB090510",
        "gbm_trigger": "bn090510016",
        "year": "2009",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": 263607749.0,
        "manual_time_window": (-5.0, 15.0),
    },
    {
        "name": "GRB130427A",
        "gbm_trigger": "bn130427324",
        "year": "2013",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB160509A",
        "gbm_trigger": "bn160509374",
        "year": "2016",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB160625B",
        "gbm_trigger": "bn160625945",
        "year": "2016",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB190114C",
        "gbm_trigger": "bn190114873",
        "year": "2019",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB090902B",
        "gbm_trigger": "bn090902462",
        "year": "2009",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB090926A",
        "gbm_trigger": "bn090926181",
        "year": "2009",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
    {
        "name": "GRB080916C",
        "gbm_trigger": "bn080916009",
        "year": "2008",
        "tte_detectors": ["n0", "n1", "n2", "n3", "n4",
                          "n5", "n6", "n7", "n8", "n9", "na", "nb"],
        "manual_trigger_met": None,
        "manual_time_window": None,
    },
]

# Base URL for Fermi/GBM TTE data at the FSSC
FSSC_BASE = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers"

# Gaussian prior on β₂ (e.g. from GW constraints or theoretical expectations)
BETA2_PRIOR_MEAN = 0.0
BETA2_PRIOR_SIGMA = 1e-4

# -------------------------------------------------------------------------
# Reduced-mode parameters
# -------------------------------------------------------------------------
# To keep the likelihood evaluation computationally manageable, we optionally
# subsample photons when there are very many in the analysis window.
MAX_PHOTONS_LIKELIHOOD = 300_000   # maximum number of photons used in likelihood
MAX_PHOTONS_PLOT = 50_000          # maximum number of photons shown in scatter plots

# -------------------------------------------------------------------------
# 1) Download TTE files for a given event
# -------------------------------------------------------------------------
def ensure_tte_files(event):
    """
    Ensure that all requested TTE files for a given GRB are present locally.

    Parameters
    ----------
    event : dict
        One entry from EVENT_CATALOG.

    Returns
    -------
    list of str
        List of local file paths to the successfully downloaded (or already
        existing) TTE files.
    """
    name = event["name"]
    trig = event["gbm_trigger"]
    year = event["year"]
    dets = event["tte_detectors"]

    out_dir = os.path.join("TTE_data", name)
    os.makedirs(out_dir, exist_ok=True)

    tte_paths = []
    for det in dets:
        fname = f"glg_tte_{det}_{trig}_v00.fit"
        url = f"{FSSC_BASE}/{year}/{trig}/current/{fname}"
        local_path = os.path.join(out_dir, fname)

        if not os.path.exists(local_path):
            print(f"  [DOWNLOAD] {name}: {url}")
            try:
                urllib.request.urlretrieve(url, local_path)
                print(f"      -> downloaded: {local_path}")
            except Exception as e:
                print(f"      -> download FAILED ({e}), skipping detector {det}")
                continue
        else:
            print(f"  [OK] {name}: already present: {local_path}")

        tte_paths.append(local_path)

    return tte_paths

# -------------------------------------------------------------------------
# 2) Mini MET-seeker: estimate trigger MET and time window
# -------------------------------------------------------------------------
def estimate_trigger_and_window(tte_files, fallback_trigger=None, fallback_window=None):
    """
    Estimate a trigger time (MET) and an analysis time window from TTE data.

    If both fallback_trigger and fallback_window are provided, they are used
    directly. Otherwise, we:
      - read all photon times from the TTE files,
      - find the minimum and maximum MET,
      - define a simple time window around that range.

    Parameters
    ----------
    tte_files : list of str
        Paths to TTE FITS files.
    fallback_trigger : float or None
        If not None, use this as the trigger MET.
    fallback_window : tuple or None
        If not None, use this as the (t_min, t_max) window relative to trigger.

    Returns
    -------
    trigger_met : float or None
        Estimated (or provided) trigger MET.
    time_window : tuple or None
        (t_min, t_max) in seconds relative to trigger_met.
    """
    if fallback_trigger is not None and fallback_window is not None:
        return fallback_trigger, fallback_window

    all_times = []
    for path in tte_files:
        try:
            with fits.open(path) as hdul:
                data = hdul["EVENTS"].data
                all_times.append(data["TIME"])
        except Exception as e:
            print(f"    [WARN] MET-seeker could not read {path} ({e})")

    if not all_times:
        print("    [ERROR] MET-seeker: no TIME data available.")
        return None, None

    times = np.concatenate(all_times)
    tmin = float(np.min(times))
    tmax = float(np.max(times))
    print(f"    [MET-seeker] MIN TIME = {tmin}, MAX TIME = {tmax}")

    # A simple heuristic: define trigger slightly before the first photon,
    # and choose a window that covers the burst duration plus some margin.
    trigger_met = tmin - 2.0
    duration = tmax - tmin
    if duration < 50.0:
        time_window = (-5.0, 40.0)
    else:
        time_window = (-5.0, min(100.0, duration + 10.0))

    print(f"    [MET-seeker] trigger_met ~ {trigger_met}, time_window ~ {time_window}")
    return trigger_met, time_window

# -------------------------------------------------------------------------
# 3) Load photons from TTE files within the chosen time window
# -------------------------------------------------------------------------
def load_photons_from_tte(tte_files, trigger_met, time_window):
    """
    Load photon arrival times and energies from TTE files, restricted to a
    given time window relative to the trigger.

    Parameters
    ----------
    tte_files : list of str
        Paths to TTE FITS files.
    trigger_met : float
        Trigger time in MET.
    time_window : tuple
        (t_min, t_max) in seconds relative to trigger_met.

    Returns
    -------
    times_rel : np.ndarray
        Photon arrival times relative to trigger (seconds).
    energies : np.ndarray
        Photon energies (instrument units or PHA).
    """
    all_times = []
    all_energies = []

    tmin_rel, tmax_rel = time_window

    for path in tte_files:
        print(f"  Loading TTE: {path}")
        try:
            with fits.open(path) as hdul:
                data = hdul["EVENTS"].data
                times = data["TIME"]

                # Use ENERGY column if present; otherwise fall back to PHA.
                if "ENERGY" in data.columns.names:
                    energies = data["ENERGY"]
                elif "PHA" in data.columns.names:
                    energies = data["PHA"].astype(float)
                else:
                    print(f"    [WARN] No ENERGY or PHA column in {path}, skipping.")
                    continue

                times_rel = times - trigger_met
                mask = (times_rel >= tmin_rel) & (times_rel <= tmax_rel)
                if np.any(mask):
                    all_times.append(times_rel[mask])
                    all_energies.append(energies[mask])
        except Exception as e:
            print(f"    [WARN] Could not read {path} ({e})")

    if not all_times:
        print("  [WARN] No photons found in the selected time window.")
        return np.array([]), np.array([])

    times_rel = np.concatenate(all_times)
    energies = np.concatenate(all_energies)
    print(f"  Total photons in window: {len(times_rel)}")
    return times_rel, energies

# -------------------------------------------------------------------------
# 4) PQG model and likelihood (with reduced mode)
# -------------------------------------------------------------------------
def pqg_time_delay(beta2, energies, ref_energy):
    """
    Quadratic PQG time delay model.

    Parameters
    ----------
    beta2 : float
        Quadratic PQG time-delay parameter β₂.
    energies : np.ndarray
        Photon energies (instrument units).
    ref_energy : float
        Reference energy (e.g. median of the sample).

    Returns
    -------
    np.ndarray
        Time delays Δt = β₂ (E - E_ref)².
    """
    return beta2 * (energies - ref_energy) ** 2

def model_emission_times(beta2, times_rel, energies, ref_energy):
    """
    Map observed arrival times to inferred emission times under the PQG model.

    t_emit = t_obs - Δt_PQG(β₂, E).

    Parameters
    ----------
    beta2 : float
        Quadratic PQG parameter.
    times_rel : np.ndarray
        Observed arrival times relative to trigger.
    energies : np.ndarray
        Photon energies.
    ref_energy : float
        Reference energy.

    Returns
    -------
    np.ndarray
        Inferred emission times.
    """
    dt = pqg_time_delay(beta2, energies, ref_energy)
    return times_rel - dt

def log_likelihood_beta2(beta2, times_rel, energies, ref_energy):
    """
    Simple likelihood for β₂ based on the variance of inferred emission times.

    The idea is that, if the PQG correction is "correct", the intrinsic
    emission times should be more tightly clustered (smaller variance).
    We model this with a Gaussian likelihood on the variance scale.

    Parameters
    ----------
    beta2 : float
        Quadratic PQG parameter.
    times_rel : np.ndarray
        Observed arrival times relative to trigger.
    energies : np.ndarray
        Photon energies.
    ref_energy : float
        Reference energy.

    Returns
    -------
    float
        Log-likelihood log p(data | β₂).
    """
    if len(times_rel) == 0:
        return -np.inf

    t_emit = model_emission_times(beta2, times_rel, energies, ref_energy)
    t_centered = t_emit - np.mean(t_emit)
    var = np.var(t_centered)

    # Effective model scale for the variance; this is a phenomenological choice.
    sigma_model = 0.1
    return -0.5 * var / (sigma_model ** 2)

def log_prior_beta2(beta2, mean, sigma):
    """
    Gaussian prior on β₂.

    Parameters
    ----------
    beta2 : float
        Quadratic PQG parameter.
    mean : float
        Prior mean.
    sigma : float
        Prior standard deviation.

    Returns
    -------
    float
        Log-prior log p(β₂).
    """
    return -0.5 * ((beta2 - mean) ** 2) / (sigma ** 2)

def compute_posterior_grid_adaptive(times_rel_full,
                                    energies_full,
                                    beta2_prior_mean,
                                    beta2_prior_sigma,
                                    beta2_min_init=-5e-4,
                                    beta2_max_init=5e-4,
                                    n_grid=1001,
                                    max_adapt=3):
    """
    Compute the posterior p(β₂ | data) on an adaptively refined 1D grid.

    Strategy:
      - Start from an initial β₂ range [beta2_min_init, beta2_max_init].
      - Evaluate the posterior on a uniform grid.
      - If the posterior maximum lies too close to a boundary, expand the range
        and repeat (up to max_adapt times).
      - Use a reduced subset of photons for the likelihood if there are more
        than MAX_PHOTONS_LIKELIHOOD photons (to keep runtime manageable).

    Parameters
    ----------
    times_rel_full : np.ndarray
        All photon arrival times relative to trigger.
    energies_full : np.ndarray
        All photon energies.
    beta2_prior_mean : float
        Prior mean for β₂.
    beta2_prior_sigma : float
        Prior standard deviation for β₂.
    beta2_min_init : float
        Initial minimum of β₂ grid.
    beta2_max_init : float
        Initial maximum of β₂ grid.
    n_grid : int
        Number of grid points.
    max_adapt : int
        Maximum number of adaptive expansions.

    Returns
    -------
    beta2_grid : np.ndarray
        Grid of β₂ values.
    posterior : np.ndarray or None
        Normalized posterior values on the grid, or None if invalid.
    mean : float
        Posterior mean of β₂.
    sigma : float
        Posterior standard deviation of β₂.
    """
    N = len(times_rel_full)
    if N == 0:
        return None, None, math.nan, math.nan

    # Reduced mode for likelihood evaluation: subsample if too many photons.
    if N > MAX_PHOTONS_LIKELIHOOD:
        print(f"  [REDUCED] Likelihood: {N} -> {MAX_PHOTONS_LIKELIHOOD} photons")
        idx = np.random.choice(N, size=MAX_PHOTONS_LIKELIHOOD, replace=False)
        times_rel = times_rel_full[idx]
        energies = energies_full[idx]
    else:
        times_rel = times_rel_full
        energies = energies_full

    ref_energy = np.median(energies)

    beta2_min = beta2_min_init
    beta2_max = beta2_max_init

    for adapt_step in range(max_adapt):
        print(
            f"  [ADAPT] β₂-grid step {adapt_step + 1}/{max_adapt}, "
            f"range=({beta2_min:.2e}, {beta2_max:.2e})"
        )

        beta2_grid = np.linspace(beta2_min, beta2_max, n_grid)
        log_post = []

        for b in tqdm(beta2_grid, desc="    β₂ grid", leave=False):
            lp = log_prior_beta2(b, beta2_prior_mean, beta2_prior_sigma)
            ll = log_likelihood_beta2(b, times_rel, energies, ref_energy)
            log_post.append(lp + ll)

        log_post = np.array(log_post)
        if not np.isfinite(log_post).any():
            print("  [ADAPT] log_posterior is -inf everywhere, aborting.")
            return beta2_grid, None, math.nan, math.nan

        # Normalize posterior in a numerically stable way.
        log_post -= np.max(log_post)
        post = np.exp(log_post)
        norm = np.trapezoid(post, beta2_grid)
        if norm == 0 or not np.isfinite(norm):
            print("  [ADAPT] Posterior normalization failed (0 or NaN), aborting.")
            return beta2_grid, None, math.nan, math.nan
        post /= norm

        mean = np.trapezoid(beta2_grid * post, beta2_grid)
        var = np.trapezoid((beta2_grid - mean) ** 2 * post, beta2_grid)
        sigma = math.sqrt(var)

        # Check whether the posterior maximum lies too close to a boundary.
        idx_max = np.argmax(post)
        frac_pos = idx_max / (len(beta2_grid) - 1)

        edge_tol = 0.02  # 2% from the edge
        if frac_pos < edge_tol:
            # Maximum at the left edge: extend range to the left.
            span = beta2_max - beta2_min
            beta2_min -= span
            print("  [ADAPT] Maximum near left edge, extending range to the left.")
            continue
        elif frac_pos > 1.0 - edge_tol:
            # Maximum at the right edge: extend range to the right.
            span = beta2_max - beta2_min
            beta2_max += span
            print("  [ADAPT] Maximum near right edge, extending range to the right.")
            continue
        else:
            # Maximum is well inside the grid: accept this range.
            print("  [ADAPT] Maximum is interior, accepting current grid.")
            return beta2_grid, post, mean, sigma

    print("  [ADAPT] Reached maximum number of adaptation steps; using last grid.")
    return beta2_grid, post, mean, sigma

# -------------------------------------------------------------------------
# 5) Plotting utilities (with subsampling for large N)
# -------------------------------------------------------------------------
def plot_beta2_posterior(beta2_grid, posterior, mean, sigma, out_png, event_name):
    """
    Plot the posterior distribution p(β₂ | data) for a single event.

    Parameters
    ----------
    beta2_grid : np.ndarray
        Grid of β₂ values.
    posterior : np.ndarray or None
        Posterior values on the grid.
    mean : float
        Posterior mean of β₂.
    sigma : float
        Posterior standard deviation of β₂.
    out_png : str
        Output path for the PNG file.
    event_name : str
        Name of the GRB event.
    """
    if posterior is None:
        print("  [WARN] No posterior available, skipping β₂ plot.")
        return

    plt.figure(figsize=(6, 4))
    plt.plot(beta2_grid, posterior, label=r"Posterior $p(\beta_2 \mid \mathrm{data})$")
    plt.axvline(
        mean,
        color="C1",
        linestyle="--",
        label=f"Mean $\\beta_2 = {mean:.2e}$",
    )
    plt.fill_between(
        beta2_grid,
        0,
        posterior,
        where=(beta2_grid >= mean - sigma) & (beta2_grid <= mean + sigma),
        color="C1",
        alpha=0.3,
        label=f"$1\\sigma$ interval: {sigma:.2e}",
    )
    plt.xlabel(r"$\beta_2$  (units consistent with PQG model)")
    plt.ylabel(r"$p(\beta_2 \mid \mathrm{GW}, \mathrm{GRB})$")
    plt.title(f"PQG quadratic time-delay parameter for {event_name}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"  Saved posterior plot to {out_png}")

def plot_toa_vs_energy(times_rel_full, energies_full, beta2_best, out_png, event_name):
    """
    Plot photon arrival times vs energy, with and without PQG correction.

    Parameters
    ----------
    times_rel_full : np.ndarray
        All photon arrival times relative to trigger.
    energies_full : np.ndarray
        All photon energies.
    beta2_best : float
        Best-fit (e.g. posterior mean) β₂ value.
    out_png : str
        Output path for the PNG file.
    event_name : str
        Name of the GRB event.
    """
    N = len(times_rel_full)
    if N == 0:
        print("  [WARN] No photons available, skipping TOA vs energy plot.")
        return

    # Subsample for plotting if there are too many photons.
    if N > MAX_PHOTONS_PLOT:
        print(f"  [REDUCED] Plot: {N} -> {MAX_PHOTONS_PLOT} points")
        idx = np.random.choice(N, size=MAX_PHOTONS_PLOT, replace=False)
        times_rel = times_rel_full[idx]
        energies = energies_full[idx]
    else:
        times_rel = times_rel_full
        energies = energies_full

    ref_energy = np.median(energies)
    t_emit = model_emission_times(beta2_best, times_rel, energies, ref_energy)

    plt.figure(figsize=(7, 5))
    plt.scatter(energies, times_rel, s=1, alpha=0.3, label="Observed $t_{\\mathrm{obs}}$")
    plt.scatter(energies, t_emit, s=1, alpha=0.3, label="PQG-corrected $t_{\\mathrm{emit}}$")
    plt.xlabel("Photon energy (instrument units)")
    plt.ylabel("Time relative to trigger (s)")
    plt.title(rf"{event_name}: arrival times vs energy with PQG correction")
    plt.legend(markerscale=5)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"  Saved TOA vs energy plot to {out_png}")

def write_summary(filename, event_name, tte_files, trigger_met, time_window,
                  beta2_prior_mean, beta2_prior_sigma,
                  n_photons, beta2_post_mean, beta2_post_sigma):
    """
    Write a human-readable text summary for a single GRB event.

    Parameters
    ----------
    filename : str
        Output path for the summary text file.
    event_name : str
        Name of the GRB event.
    tte_files : list of str
        List of TTE file paths used.
    trigger_met : float
        Trigger time in MET.
    time_window : tuple
        (t_min, t_max) in seconds relative to trigger.
    beta2_prior_mean : float
        Prior mean for β₂.
    beta2_prior_sigma : float
        Prior standard deviation for β₂.
    n_photons : int
        Number of photons used in the analysis window.
    beta2_post_mean : float
        Posterior mean of β₂.
    beta2_post_sigma : float
        Posterior standard deviation of β₂.
    """
    with open(filename, "w") as f:
        f.write("PQG vs GRB photon-channel pipeline summary\n")
        f.write("=========================================\n\n")
        f.write(f"Generated at: {datetime.now(UTC).isoformat()} UTC\n\n")

        f.write("Input configuration:\n")
        f.write(f"  Event name: {event_name}\n")
        f.write(f"  TTE files: {', '.join(os.path.basename(p) for p in tte_files)}\n")
        f.write(f"  Trigger MET: {trigger_met}\n")
        f.write(
            f"  Time window: [{time_window[0]}, {time_window[1]}] s "
            "relative to trigger\n"
        )
        f.write(
            f"  Prior on beta2 (from GW or theory): "
            f"mean={beta2_prior_mean}, sigma={beta2_prior_sigma}\n\n"
        )

        f.write("Results (PQG photon-channel analysis):\n")
        f.write(f"  Number of photons used: {n_photons}\n")
        f.write("  Used energy information: True\n")
        f.write(f"  Posterior beta2 mean: {beta2_post_mean}\n")
        f.write(f"  Posterior beta2 sigma: {beta2_post_sigma}\n\n")

        f.write("Interpretation (to be expanded in the paper):\n")
        f.write("  - Here one can discuss how this constraint compares to the\n")
        f.write("    GW-only (LIGO/Virgo/KAGRA) constraints and to other channels,\n")
        f.write("    and what it implies for Planck-scale quantum gravity models.\n")

    print(f"  Summary written to {filename}")

# -------------------------------------------------------------------------
# 6) Run the analysis for a single event
# -------------------------------------------------------------------------
def run_event(event_cfg):
    """
    Run the full PQG photon-channel analysis for a single GRB event.

    Steps:
      - ensure TTE files are available,
      - estimate trigger MET and time window (or use manual values),
      - load photons in that window,
      - compute the β₂ posterior on an adaptive grid,
      - generate plots and a text summary.

    Parameters
    ----------
    event_cfg : dict
        One entry from EVENT_CATALOG.

    Returns
    -------
    dict or None
        Dictionary with event name, number of photons, β₂ mean and σ, or
        None if the event could not be analyzed.
    """
    name = event_cfg["name"]
    print(f"\n=== Running PQG photon-channel analysis for {name} ===")

    tte_files = ensure_tte_files(event_cfg)
    if not tte_files:
        print("  [SKIP] No TTE files available for this event.")
        return None

    trig_manual = event_cfg.get("manual_trigger_met", None)
    win_manual = event_cfg.get("manual_time_window", None)
    trigger_met, time_window = estimate_trigger_and_window(
        tte_files, fallback_trigger=trig_manual, fallback_window=win_manual
    )
    if trigger_met is None or time_window is None:
        print("  [SKIP] Could not determine trigger_met/time_window.")
        return None

    times_rel_full, energies_full = load_photons_from_tte(
        tte_files, trigger_met, time_window
    )
    if len(times_rel_full) == 0:
        print("  [SKIP] 0 photons in the analysis window.")
        return None

    beta2_grid, post, mean, sigma = compute_posterior_grid_adaptive(
        times_rel_full,
        energies_full,
        BETA2_PRIOR_MEAN,
        BETA2_PRIOR_SIGMA,
        beta2_min_init=-5e-4,
        beta2_max_init=5e-4,
        n_grid=1001,
        max_adapt=3,
    )
    if post is None or not np.isfinite(sigma):
        print("  [SKIP] Posterior is not well-defined for this event.")
        return None

    out_dir = os.path.join(BASE_OUTPUT_DIR, name)
    os.makedirs(out_dir, exist_ok=True)

    posterior_png = os.path.join(out_dir, f"beta2_posterior_{name}.png")
    toa_png = os.path.join(out_dir, f"toa_vs_energy_{name}.png")
    summary_txt = os.path.join(out_dir, f"summary_{name}_PQG.txt")

    plot_beta2_posterior(beta2_grid, post, mean, sigma, posterior_png, name)
    plot_toa_vs_energy(times_rel_full, energies_full, mean, toa_png, name)
    write_summary(
        summary_txt,
        name,
        tte_files,
        trigger_met,
        time_window,
        BETA2_PRIOR_MEAN,
        BETA2_PRIOR_SIGMA,
        n_photons=len(times_rel_full),
        beta2_post_mean=mean,
        beta2_post_sigma=sigma,
    )

    return {
        "name": name,
        "n_photons": len(times_rel_full),
        "beta2_mean": mean,
        "beta2_sigma": sigma,
    }

# -------------------------------------------------------------------------
# 7) Global combination of β₂ constraints across events
# -------------------------------------------------------------------------
def combine_beta2(results):
    """
    Combine β₂ constraints from multiple events using inverse-variance weighting.

    Parameters
    ----------
    results : list of dict
        Each dict must contain 'beta2_mean' and 'beta2_sigma'.

    Returns
    -------
    beta2_combined : float
        Inverse-variance–weighted mean of β₂.
    sigma_combined : float
        Corresponding combined standard deviation.
    """
    weights = []
    means = []
    for r in results:
        if r["beta2_sigma"] > 0 and np.isfinite(r["beta2_sigma"]):
            w = 1.0 / (r["beta2_sigma"] ** 2)
            weights.append(w)
            means.append(r["beta2_mean"])

    if not weights:
        return math.nan, math.nan

    weights = np.array(weights)
    means = np.array(means)
    beta2_combined = np.sum(weights * means) / np.sum(weights)
    sigma_combined = math.sqrt(1.0 / np.sum(weights))
    return beta2_combined, sigma_combined

# -------------------------------------------------------------------------
# 8) Main driver
# -------------------------------------------------------------------------
def main():
    """
    Main entry point for the PQG photon-channel multi-event pipeline.

    Loops over all events in EVENT_CATALOG, runs the analysis for each,
    writes a global CSV table, and prints the combined β₂ constraint and
    its Z-score.
    """
    all_results = []

    for ev in EVENT_CATALOG:
        try:
            res = run_event(ev)
            if res is not None:
                all_results.append(res)
        except Exception as e:
            print(f"  [ERROR] Event {ev['name']} failed with exception: {e}")
            continue

    # Write global results table
    csv_path = os.path.join(BASE_OUTPUT_DIR, "beta2_results_all_events.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["event_name", "n_photons", "beta2_mean", "beta2_sigma"])
        for r in all_results:
            writer.writerow(
                [r["name"], r["n_photons"], r["beta2_mean"], r["beta2_sigma"]]
            )
    print(f"\nGlobal results table written to {csv_path}")

    # Combine β₂ constraints across all events
    beta2_comb, sigma_comb = combine_beta2(all_results)
    print("\n=== Combined beta2 over photon-channel events (inverse-variance weighted) ===")
    print(f"  beta2_combined = {beta2_comb:.3e}")
    print(f"  sigma_combined = {sigma_comb:.3e}")
    if np.isfinite(beta2_comb) and np.isfinite(sigma_comb) and sigma_comb > 0:
        Z = beta2_comb / sigma_comb
        print(f"  Z = {Z:.2f} sigma")
    else:
        print("  Z is not well-defined (insufficient valid events).")

if __name__ == "__main__":
    main()
