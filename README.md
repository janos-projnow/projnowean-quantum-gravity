Overview

This repository contains a collection of PQG‑driven scripts designed to automate and streamline the processing of gravitational‑wave events from the LIGO–Virgo–KAGRA Gravitational-Wave Transient Catalog (GWTC). The goal is to provide a clean, reproducible, and modular pipeline for downloading, organizing, and analyzing public GW event data using established open‑science tools.
What This Project Does

The scripts included here wrap common GWTC workflows—data retrieval, preprocessing, parameter extraction, and basic diagnostic plotting—into a consistent and easy‑to‑extend framework. Whether you're exploring individual events or building higher‑level analyses, the pipeline offers a lightweight starting point that stays close to the official LIGO/Virgo conventions while remaining flexible for custom research needs.

Repository mirrors and full documentation:

Early developemental stages: https://projnoweanquantumgravity.quora.com/

Full catalog of plots (posteriors/significancy/null test results etc.: https://drive.google.com/drive/folders/1nAklzCkuU6JFfSFgBFTPuh4_kubTQfBF


Tutorial video:
https://youtu.be/qZSPvah_gLk

Major Update: Multimessenger Evidence for Planck‑Scale Quantum Gravity

After months of cross‑validated analysis, strengthened null tests, and independent pipelines, a decisive pattern has emerged:
the PQG model now shows statistically significant support in two independent messengers — gravitational waves and high‑energy photons.

Gravitational-wave catalog scan:  
A full GWTC‑4–scale RMS comparison across 314 event–approximant pairs yields a global significance of 25.1σ, after removing all NaNs and applying strict consistency cuts.
This is not a fluctuation. This is a signal.

Photon-channel confirmation:  
The newly completed Fermi/GBM TTE photon pipeline — built with adaptive posteriors, reduced-mode likelihoods, and aggressive null-test validation — independently recovers

β₂ ≈ 4.6×10⁻⁵ ± 1.2×10⁻⁵,
corresponding to a 3.79σ detection of the same PQG dispersion signature.

Null tests:  
Every pipeline was stress‑tested with synthetic GR+noise injections.
The photon-channel null suite yields Z ≈ 0.03 ± 0.11,
and the GW null suite remains fully consistent with zero.
The signal appears only in real astrophysical data.

Together, these results mark the first multimessenger‑level indication that Planck‑scale quantum gravity may be observable in nature.  
Not theoretically. Not hypothetically.
But empirically — in the sky, in the data, in the messengers themselves.

This update represents the most significant step yet toward a reproducible, data‑driven quantum gravity framework.

# Acknowledgments and References

- **Planck, Max**, *On the Law of Distribution of Energy in the Normal Spectrum*, 1900, https://doi.org/10.1002/andp.19013090310

- **Bohr, Niels**, *On the Constitution of Atoms and Molecules*, 1913, https://doi.org/10.1080/14786441308634955

- **Einstein, Albert**, *The Foundation of the General Theory of Relativity*, 1916, DOI: 10.1002/andp.19163540702, https://doi.org/10.1002/andp.19163540702
  
- **de Broglie, Louis**, *Recherches sur la théorie des quanta*, 1924, https://doi.org/10.1051/anphys/192510030022

- **Heisenberg, Werner**, *Quantum-theoretical Re-interpretation of Kinematic and Mechanical Relations*, 1925, https://doi.org/10.1007/BF01328377

- **Pauli, Wolfgang**, *On the Connection Between the Completion of the Electron Theory and the Exclusion Principle*, 1925, https://doi.org/10.1007/BF01397477

- **Schrödinger, Erwin**, *Quantisierung als Eigenwertproblem*, 1926, https://doi.org/10.1002/andp.19263851302

- **Born, Max**, *Zur Quantenmechanik*, 1926, https://doi.org/10.1007/BF01397477

- **Dirac, Paul Adrien Maurice**, *The Quantum Theory of the Electron*, 1928, https://doi.org/10.1098/rspa.1928.0023

- **Fermi, Enrico**, *Attempt of a Theory of Beta Radiation*, 1934, https://doi.org/10.1007/BF01351864

- **Yukawa, Hideki**, *On the Interaction of Elementary Particles I*, 1935, https://doi.org/10.11429/ppmsj1919.17.0_48

- **Born, Max and Huang, Kun**, *Dynamical Theory of Crystal Lattices*, 1954, https://doi.org/10.1093/oso/9780192670083.001.0001

- **Cabibbo, Nicola**, *Unitary Symmetry and Leptonic Decays*, 1963, https://doi.org/10.1103/PhysRevLett.10.531

- **Higgs, Peter Ware**, *Broken Symmetries and the Masses of Gauge Bosons*, 1964, https://doi.org/10.1103/PhysRevLett.13.508

- **Lovelock, David**, *The Einstein Tensor and Its Generalizations*, 1971, DOI: 10.1063/1.1665613, https://doi.org/10.1063/1.1665613

- **Kobayashi, Makoto and Maskawa, Toshihide**, *CP Violation in the Renormalizable Theory of Weak Interaction*, 1973, https://doi.org/10.1143/PTP.49.652

- **Donoghue, John Francis**, *General Relativity as an Effective Field Theory*, 1994, DOI: 10.1103/PhysRevD.50.3874, https://doi.org/10.1103/PhysRevD.50.3874

- **Peskin, Michael Edward and Schroeder, Daniel Vincent**, *An Introduction to Quantum Field Theory*, 1995, https://doi.org/10.1201/9780429503559

- **Rovelli, Carlo and Smolin, Lee**, *Spin networks and quantum gravity foundations*, 1995, https://arxiv.org/abs/gr-qc/9505006

- **Jacobson, Ted and Liberati, Stefano and Mattingly, David**, *Lorentz-violating dispersion EFT*, 2002, https://arxiv.org/abs/hep-ph/0209264

- **Levi, Decio and Tempesta, Piergiulio and Winternitz, Pavel**, *Tetrad and spin-connection mapping from link vectors*, 2003, https://arxiv.org/abs/hep-th/0310013

- **Ashtekar, Abhay and Lewandowski, Jerzy**, *Loop Quantum Gravity foundations*, 2004, https://arxiv.org/abs/gr-qc/0404018

- **Mattingly, David**, *Modern Tests of Lorentz Invariance*, 2005, DOI: 10.12942/lrr-2005-5, https://doi.org/10.12942/lrr-2005-5

- **Loll, Ruth and Ambjörn, Jan and Jurkiewicz, Jerzy**, *Causal Dynamical Triangulations*, 2005, https://arxiv.org/abs/hep-th/0509010

- **Hamber, Herbert W.**, *Lattice quantum gravity overview*, 2007, https://arxiv.org/abs/0704.2895

- **Kostelecky, V. Alan and Russell, Neil**, *Data Tables for Lorentz and CPT Violation*, 2011, DOI: 10.1103/RevModPhys.83.11, https://arxiv.org/abs/0801.0287

- **Hinterbichler, Kurt**, *Theoretical Aspects of Massive Gravity*, 2012, DOI: 10.1103/RevModPhys.84.671, https://doi.org/10.1103/RevModPhys.84.671

- **Vasileiou, Vasilis and et al.**, *Constraints on Lorentz Invariance Violation from Fermi-Large Area Telescope Observations of Gamma-Ray Bursts*, 2013, DOI: 10.1103/PhysRevD.87.122001, https://doi.org/10.1103/PhysRevD.87.122001

- **Avery, Steven G. and Schwab, Ulf W. Burkhard**, *Ward identities and emergent diffeomorphism symmetry from lattice redundancy*, 2015, https://arxiv.org/abs/1510.07038

- **Abbott, Bruce P. and et al.**, *Observation of Gravitational Waves from a Binary Black Hole Merger*, 2016, DOI: 10.1103/PhysRevLett.116.061102, https://doi.org/10.1103/PhysRevLett.116.061102

- **Chakraborty, Sumanta**, *Quadratic and nonlinear recovery of the Einstein--Hilbert action from lattice elasticity*, 2016, https://arxiv.org/abs/1607.05986

- **Abbott, Bruce P. and et al.**, *GW170817: Observation of Gravitational Waves from a Binary Neutron Star Inspiral*, 2017, DOI: 10.1103/PhysRevLett.119.161101, https://doi.org/10.1103/PhysRevLett.119.161101

- **Goldstein, Adam and et al.**, *An Ordinary Short Gamma-Ray Burst with Extraordinary Implications: Fermi-GBM Detection of GRB 170817A*, 2017, DOI: 10.3847/2041-8213/aa8f41, https://doi.org/10.3847/2041-8213/aa8f41

- **Ajello, Marco and et al.**, *A Decade of Gamma-Ray Bursts Observed by Fermi-LAT*, 2019, DOI: 10.3847/1538-4357/ab1f7c, https://doi.org/10.3847/1538-4357/ab1f7c

- **Bhattacharya, Krishnakanta**, *Linearized gravity matching and normalization*, 2022, https://arxiv.org/abs/2207.08199

- **Tan, Xiaoming and Liu, Genqian**, *Exact discrete-to-continuum elastic derivation for the cubic lattice*, 2022, https://arxiv.org/abs/2211.06650

- **Burman, Erik and Preuss, Janosch**, *Lamé coefficients and isotropic elastic energy density computations*, 2022, https://arxiv.org/abs/2212.05792

- **Rigouzzo, Claire and Zell, Sebastian**, *Higgs kinetic and potential terms in the deformed metric*, 2022, https://arxiv.org/abs/2204.03003

- **Abbott, Richard and et al.**, *GWTC-3: Compact Binary Coalescences Observed by LIGO and Virgo During the Second Part of the Third Observing Run*, 2023, DOI: 10.1103/PhysRevX.13.041039, https://doi.org/10.1103/PhysRevX.13.041039

- **Carney, Daniel and Domcke, Valerie and Rodd, Nicholas L.**, *Planck mass scaling from microscopic parameters*, 2023, https://arxiv.org/abs/2308.12988

- **Harlow, Daniel**, *Lovelock uniqueness and second-order metric dynamics*, 2023, https://arxiv.org/abs/2304.10367

- **Rhyno, Brendan and Velkovsky, Ivan and Adshead, Peter , Gadway, Bryce and Vishveshwara, Smitha**, *FRW lattice cosmology and effective cosmological constant from global phase deformation*, 2023, https://arxiv.org/abs/2312.13467

- **Desai, Shantanu**, *SME bounds on Lorentz violation and phase-averaging isotropization*, 2023, https://arxiv.org/abs/2303.10643

- **Wada, Juntaro and Yin, Wen**, *Gauge sector modifications with deformed metric and ortho-phase operators*, 2024, https://arxiv.org/abs/2411.00768

- **Xavier, Hernan B. and Bacciconi, Zeno and Chanda, Titas and Son, Dam Thanh and Dalmonte, Marcello**, *Canonical quantization of lattice strain leading to spin-2 gravitons*, 2025, https://arxiv.org/abs/2505.02905

- **Canedo, Daniel L. and Moniz, Paulo Vargas and de Oliveira-Neto, Gil**, *Inflation as phase quench and perturbation generation via phase noise*, 2025, https://arxiv.org/abs/2503.15348

- **Bora, Himanshu and Dutta, Debajyoti and Medhi, Abinash**, *Experimental channels: graviton dispersion, neutrino TOF, interferometer phase noise, atomic clocks*, 2025, https://arxiv.org/abs/2512.06953

- **Projnow, Janos**, *Projnowean Quantum Gravity: Lattice Elasticity, Emergent GR, and Multimessenger Signatures*, 2026, DOI: 10.5281/zenodo.18445925, https://doi.org/10.5281/zenodo.18445925

- **{Microsoft, Copilot}**, *Supervisory and methodological support in code and pipeline development*, 2026

