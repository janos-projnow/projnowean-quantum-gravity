================================================================================
          PROJNOWEAN QUANTUM GRAVITY (PQG) REPRODUCIBILITY SUITE
================================================================================

This repository

https://drive.google.com/drive/folders/1JBdwQJhQootK-xKx2_FYqfRw-3Stex0W

contains the official Singularity / Apptainer image (.sif) and 
pipeline scripts for the Projnowean Quantum Gravity (PQG) data analysis, 
including gravitational wave, gamma-ray (Fermi/LAT photon timing), and null-test 
validations.

--------------------------------------------------------------------------------
1. OVERVIEW & INCLUDED ANALYSIS SCRIPTS
--------------------------------------------------------------------------------
The pipeline consists of 4 main Python scripts pre-installed inside the image:

  1) PQG_photon_channel.py
     Pertains to Fermi/LAT photon timing & quantum gravity dispersion analysis.
  2) PQG_FULL_CATALOG_SCAN_PLOT.py
     Executes full catalog gravitational wave (LVK) data fitting and plotting.
  3) null_test.py
     Generates statistical null-hypothesis tests for GW data.
  4) foton_null_test.py
     Generates statistical null-hypothesis tests for photon channel data.

--------------------------------------------------------------------------------
2. PREREQUISITES & COMPATIBILITY NOTES
--------------------------------------------------------------------------------
- Requires Apptainer >= 1.0.0 or Singularity >= 3.8.0 installed on the host system.
- Host Operating Systems: Compatible with Linux distributions (Ubuntu 20.04+, 
  Debian 11+, RHEL/CentOS 8+, Arch Linux, Fedora).
- Read-Only Container Filesystem: By default, the container environment is 
  read-only. Output directories (such as 'PQG_Photon_channel/') must be written 
  to the host system. This is handled via directory binding (--bind .:/opt/pqg).

--------------------------------------------------------------------------------
3. EXECUTION METHODS & ALL-IN-ONE RUNS
--------------------------------------------------------------------------------

Below are the two recommended methods to execute all 4 scripts sequentially in a 
single ("all-in-one") execution pipeline.

================================================================================
METHOD 1: STANDARD DIRECTORY BINDING (RECOMMENDED)
================================================================================
Use this method for most standard Linux installations where user namespaces and 
fuse mounts are properly configured.

Execution command (All-in-One Sequential Run):

apptainer exec --bind .:/opt/pqg --pwd /opt/pqg pqg_pipeline.sif bash -c \
  "python3 PQG_photon_channel.py && \
   python3 PQG_FULL_CATALOG_SCAN_PLOT.py && \
   python3 null_test.py && \
   python3 foton_null_test.py"

--------------------------------------------------------------------------------
METHOD 2: WRITABLE OVERLAY BINDING (FOR STRICT / RESTRICTED ENVIRONMENTS)
================================================================================
Use this method if your host environment enforces strict read-only filesystems, 
lacks full fuse3 permissions, or if scripts attempt to write temporary files 
to system locations like /tmp inside the container.

Execution command (All-in-One Sequential Run):

apptainer exec --writable-tmpfs --bind .:/opt/pqg --pwd /opt/pqg pqg_pipeline.sif bash -c \
  "python3 PQG_photon_channel.py && \
   python3 PQG_FULL_CATALOG_SCAN_PLOT.py && \
   python3 null_test.py && \
   python3 foton_null_test.py"

================================================================================
4. VERIFICATION & OUTPUTS
================================================================================
Upon successful completion of the all-in-one execution, all generated datasets, 
fit parameters, and high-resolution diagnostic plots will be stored directly in 
your working directory on the host machine.

To quick-test container environment integrity before running the full suite:

  apptainer exec pqg_pipeline.sif python3 -c "import numpy, scipy, gwpy, pycbc; print('OK')"

If 'OK' is printed, the container environment and scientific libraries are ready.
================================================================================