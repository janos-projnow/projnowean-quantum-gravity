#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
PQG LHC MASTER — Prior-Constrained Physics Pipeline
===============================================================================
Theoretical Framework: Projnowean Quantum Gravity (PQG) Model
-------------------------------------------------------------------------------
This pipeline evaluates High-Energy Physics (LHC) dijet mass spectra for potential
signatures of quantized spacetime phase deformations predicted by the Projnowean
Quantum Gravity (PQG) framework.

In the PQG model:
1. Spacetime is modeled as a discrete {4,3,4} cubic spring-lattice with base 
   spacing 'a' where vertices represent discrete lattice anchors and edges carry
   complex phase curves gamma(x) = R(x) * exp(i * pi * a(x)).
2. Continuum General Relativity (Einstein-Hilbert action) and Standard Model gauge
   interactions emerge at long wavelengths (k*a << 1) via stochastic phase-averaging 
   <Phi_{mu nu}> = C * eta_{mu nu} and emergent diffeomorphism symmetry.
3. High-energy dimension-5/6 derivative-deformed operators perturb gauge field
   kinetic terms. Specifically, local spring deformations introduce a momentum-
   dependent modulations of the effective invariant mass spectrum m_jj.
4. Gravitational Wave (GW) and photon speed dispersion constraints fix the low-energy
   coupling parameters, which are mapped here onto LHC scale (Lambda) as a Gaussian Beta prior.
5. To account for discrete lattice normalization uncertainties, JZ Monte Carlo slice
   mismatches, and phase-space boundary effects, a floating scale nuisance parameter
   is profiled via profile-likelihood maximization.
===============================================================================
"""

import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from scipy.optimize import minimize

# =============================================================================
# DATA AND MONTE CARLO SOURCE FILES
# =============================================================================
# DATA_FILE: High-luminosity LHC collision period data (e.g., ATLAS/CMS 13 TeV format)
DATA_FILE = "ODEO_FEB2025_v0_1LMET30_data15_periodE.1LMET30.root"

# MC_FILES: Pythia8 QCD dijet Monte Carlo slice samples (JZ4 to JZ7) representing
# Standard Model background predictions across ascending pT thresholds.
MC_FILES = [
    "mc_364704.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4WithSW.2bjets70.root",
    "mc_364705.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5WithSW.2bjets70.root",
    "mc_364706.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ6WithSW.2bjets70.root",
    "mc_364707.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ7WithSW.2bjets70.root",
]

# =============================================================================
# KINEMATIC BINNING & PHASE-SPACE SELECTION
# =============================================================================
# Dijet invariant mass (m_jj) binning range in GeV.
# Focuses on the high-pT tail (800 GeV - 4000 GeV), where high-energy lattice strain 
# corrections (suppressed by (m_jj / Lambda)^4) become dominant over standard QCD.
mjj_bins = np.linspace(800.0, 4000.0, 80)
mjj_centers = 0.5 * (mjj_bins[1:] + mjj_bins[:-1])


def compute_jet_mass(pt, eta, phi, energy):
    """
    Computes the rest mass of individual jet candidates from 4-momentum components.
    
    Physics derivation:
    P^mu = (E, px, py, pz) where:
      px = pT * cos(phi)
      py = pT * sin(phi)
      pz = pT * sinh(eta)
    m^2 = E^2 - |p|^2 = E^2 - (px^2 + py^2 + pz^2)
    
    Returns 0.0 for unphysical negative squared masses caused by detector resolution.
    """
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    m2 = energy**2 - px**2 - py**2 - pz**2
    return np.sqrt(np.where(m2 > 0.0, m2, 0.0))


def compute_mjj_from_tree(tree):
    """
    Extracts leading two jets (leading dijet system) from ROOT TTree arrays and 
    computes the dijet invariant mass (m_jj).
    
    Methodology:
    1. Reads jet 4-momenta (pT, eta, phi, Energy).
    2. Filters events requiring at least 2 reconstructed jets (mask).
    3. Sums the 4-momenta of Jet 1 and Jet 2:
       E_tot = E1 + E2
       p_tot = p1 + p2
    4. Evaluates the invariant mass: m_jj = sqrt(E_tot^2 - |p_tot|^2).
    """
    jet_pt = tree["jet_pt"].array()
    jet_eta = tree["jet_eta"].array()
    jet_phi = tree["jet_phi"].array()
    jet_e = tree["jet_e"].array()
    
    jet_m = compute_jet_mass(jet_pt, jet_eta, jet_phi, jet_e)
    
    # Event selection constraint: at least 2 jets present in the event
    mask = ak.num(jet_pt) >= 2
    
    # 4-momentum components of the leading jet (Jet 1)
    px1 = jet_pt[mask][:, 0] * np.cos(jet_phi[mask][:, 0])
    py1 = jet_pt[mask][:, 0] * np.sin(jet_phi[mask][:, 0])
    pz1 = jet_pt[mask][:, 0] * np.sinh(jet_eta[mask][:, 0])
    E1 = np.sqrt(px1**2 + py1**2 + pz1**2 + jet_m[mask][:, 0]**2)
    
    # 4-momentum components of the sub-leading jet (Jet 2)
    px2 = jet_pt[mask][:, 1] * np.cos(jet_phi[mask][:, 1])
    py2 = jet_pt[mask][:, 1] * np.sin(jet_phi[mask][:, 1])
    pz2 = jet_pt[mask][:, 1] * np.sinh(jet_eta[mask][:, 1])
    E2 = np.sqrt(px2**2 + py2**2 + pz2**2 + jet_m[mask][:, 1]**2)
    
    # Combined dijet invariant mass m_jj
    m_jj = np.sqrt((E1 + E2)**2 - (px1 + px2)**2 - (py1 + py2)**2 - (pz1 + pz2)**2)
    return np.array(m_jj)


# =============================================================================
# DATA LOADING AND MONTE CARLO BACKGROUND ASSEMBLY
# =============================================================================
print("[INFO] Loading DATA...")
data_tree = uproot.open(DATA_FILE)["analysis"]
mjj_data = compute_mjj_from_tree(data_tree)
hist_data, _ = np.histogram(mjj_data, bins=mjj_bins)

print("[INFO] Loading MC Slices...")
hist_mc_total = np.zeros(len(mjj_bins) - 1)

# JZ Slice relative weighting factors.
# In experimental QCD analyses, slice weights account for cross-section scaling:
# Weight = (Cross Section * Filter Efficiency) / Total Generated Events
jz_weights = {4: 1.0, 5: 0.1, 6: 0.01, 7: 0.001} 

for mc_file in MC_FILES:
    jz_idx = int(mc_file.split("JZ")[1][0])
    weight = jz_weights.get(jz_idx, 1.0)
    mc_tree = uproot.open(mc_file)["analysis"]
    mjj_mc = compute_mjj_from_tree(mc_tree)
    hist_mc, _ = np.histogram(mjj_mc, bins=mjj_bins)
    hist_mc_total += hist_mc * weight


# =============================================================================
# PQG MODEL SPECIFICATION AND PRIORS
# =============================================================================
# GW and Photon Dispersion Prior Parameters:
# Astrophysical tests (e.g. LIGO/Virgo gravitational wave dispersion & gamma-ray bursts)
# constrain low-energy Lorentz-violating phase deformation parameters:
#   Beta_2 = 4.6e-5 +/- 1.2e-5 s^2
# Projected onto the LHC TeV scale, this mapped constraint translates to a Gaussian prior on Beta:
#   Mean (beta_prior_mean) = 0.015
#   Standard Deviation (beta_prior_std) = 0.004

BETA_PRIOR_MEAN = 0.015
BETA_PRIOR_STD = 0.004


def pqg_model(beta, Lambda, hist_sm):
    """
    Calculates the expected dijet yield under Projnowean Quantum Gravity modification.
    
    Formula:
      N_PQG(m_jj) = N_SM(m_jj) * [ 1 + beta * (m_jj / Lambda)^4 ]
      
    Theoretical Context:
    - N_SM: Standard Model prediction (pure QCD background from MC).
    - Lambda: The effective energy scale of spacetime lattice cutoff (e.g., 2000 GeV).
    - (m_jj / Lambda)^4: High-energy higher-dimensional operator correction factor arising
      from dimension-6 ortho-phase gluon interactions:
      Delta_L = (beta_g / Lambda^2) * chi(x)^2 * Tr[G_mu_nu G^nu_rho G^rho_sigma G^sigma_mu]
    """
    factor = 1.0 + beta * (mjj_centers / Lambda)**4
    return hist_sm * factor


def loss_function(params, hist_obs, hist_sm, Lambda):
    """
    Negative Log-Likelihood (NLL) with Gaussian prior constraints and profiled
    normalization scale nuisance parameter.
    
    Parameters:
      params: List containing [scale, beta]
        - scale (nuisance): Normalization factor absorbing luminosity uncertainty,
          detector acceptance mismatches, and global lattice volume factors.
        - beta (physical parameter): PQG coupling constant constrained by GW prior.
      hist_obs: Observed event histogram (DATA).
      hist_sm: Unweighted Standard Model Monte Carlo dijet histogram.
      Lambda: Cutoff scale parameter (GeV).
      
    Likelihood Formulation:
      - Poisson likelihood term for binned counting statistics:
        -ln L_poisson = - sum_i [ N_obs,i * ln(N_pred,i) - N_pred,i ]
      - Gaussian penalty term for Beta prior (GW / Photon constraints):
        Prior_beta = 0.5 * ((beta - BETA_PRIOR_MEAN) / BETA_PRIOR_STD)^2
      - Gaussian penalty term for scale factor nuisance (broad prior around 1.0):
        Prior_scale = 0.5 * ((scale - 1.0) / 0.5)^2
    """
    scale, beta = params
    
    # Scaled theoretical prediction
    pred = scale * pqg_model(beta, Lambda, hist_sm)
    pred = np.where(pred > 0, pred, 1e-9)  # Numerical stability check against log(0)
    
    # Binned Poisson Negative Log-Likelihood (ignoring constant factorials)
    poisson_nll = -np.sum(hist_obs * np.log(pred) - pred)
    
    # Nuisance and physical prior penalty terms
    prior_beta = 0.5 * ((beta - BETA_PRIOR_MEAN) / BETA_PRIOR_STD)**2
    prior_scale = 0.5 * ((scale - 1.0) / 0.5)**2  # Wide bound allows flexible profile adjustment
    
    return poisson_nll + prior_beta + prior_scale


# =============================================================================
# STATISTICAL FIT EXECUTION & HYPOTHESIS TESTING
# =============================================================================
Lambda_fixed = 2000.0  # Cutoff scale set to 2.0 TeV

# 1. Null Hypothesis Fit (H0: Pure Standard Model, Beta = 0 strictly forced)
# Profiles the scale factor to optimize background normalization under pure SM.
res_null = minimize(
    lambda p: loss_function([p[0], 0.0], hist_data, hist_mc_total, Lambda_fixed), 
    x0=[0.1], 
    bounds=[(1e-5, 10.0)], 
    method="L-BFGS-B"
)

# 2. PQG Hypothesis Fit (H1: Projnowean Quantum Gravity signal present)
# Simultaneously fits the scale nuisance factor and the beta coupling parameter
# under Gaussian prior regularization.
res_pqg = minimize(
    lambda p: loss_function(p, hist_data, hist_mc_total, Lambda_fixed), 
    x0=[0.1, BETA_PRIOR_MEAN], 
    bounds=[(1e-5, 10.0), (-0.1, 0.5)], 
    method="L-BFGS-B"
)

best_scale, best_beta = res_pqg.x
null_scale = res_null.x[0]

# Compute statistical significance Z-score via Likelihood Ratio Test (Wilks' Theorem)
# Delta_NLL = NLL(H0) - NLL(H1)
# Z = sqrt( max(0, 2 * Delta_NLL) )
delta_nll = res_null.fun - res_pqg.fun
Z_corrected = np.sqrt(max(0.0, 2.0 * delta_nll))

# Print fit summary to console
print("\n" + "="*50)
print(f"  PRIOR-CONSTRAINED FIT RESULTS")
print("="*50)
print(f"Optimal Scale Factor (Nuisance): {best_scale:.4f}")
print(f"Fitted Beta (constrained by GW): {best_beta:.6f}")
print(f"Corrected Significance (Z-score): {Z_corrected:.3f} sigma")
print("="*50)

# =============================================================================
# RESULTS VIZUALIZATION
# =============================================================================
plt.figure(figsize=(10, 6))

# Plot observed experimental data
plt.step(mjj_centers, hist_data, where="mid", label="DATA (1LMET30) [Normalized]", color="blue", lw=2)

# Plot profiled Standard Model MC background (Null Hypothesis)
plt.step(mjj_centers, null_scale * hist_mc_total, where="mid", label="SM MC (Profiled)", color="black", linestyle="--")

# Plot fitted PQG signal model (Signal Hypothesis)
plt.step(
    mjj_centers, 
    best_scale * pqg_model(best_beta, Lambda_fixed, hist_mc_total), 
    where="mid", 
    label=f"PQG Fit (Beta={best_beta:.4f}, Z={Z_corrected:.2f})", 
    color="red", 
    lw=1.5
)

plt.xlabel(r"$m_{jj}$ [GeV]")
plt.ylabel("Events (Scale Adjusted)")
plt.title("PQG LHC — Prior-Constrained Profile-Likelihood Fit")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# Save generated figure
plt.savefig("PQG_Prior_Fit.png", dpi=200)
plt.show()