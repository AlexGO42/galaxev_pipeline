#!/bin/bash

# An example batch script that calls stellar_photometrics.py with an
# appropriate set of parameters.

# Any user-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# A folder to store all the output corresponding to this set of synthetic
# images (with subdirectories corresponding to different simulations,
# snapshots, etc., as will become clear later):
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Simulation snapshot, although not really needed yet (only for the filename):
SNAPNUM=78

# Observation redshift, which usually corresponds to snapnum:
USE_Z=0.297717684517

# Directory with the GALAXEV model files:
BC03_MODEL_DIR=/path/to/bc03/Padova1994/chabrier

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Illustris or IllustrisTNG (other simulations eventually):
SUITE=IllustrisTNG

# Calculate magnitudes with or without Charlot & Fall (2000) dust model:
USE_CF00=0  # 0 = no, 1 = yes

python ${CODEDIR}/galaxev_pipeline/stellar_photometrics.py \
    ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} \
    ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_SET}
