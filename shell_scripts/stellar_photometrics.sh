#!/bin/bash

# An example batch script that calls stellar_photometrics.py with an
# appropriate set of parameters.

# Any user-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# A folder to store all the output corresponding to this set of synthetic
# images (with subdirectories corresponding to different simulations,
# snapshots, etc., as will become clear later):
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Directory with the GALAXEV model files:
BC03_MODEL_DIR=/path/to/bc03/Padova1994/chabrier

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Specify simulation and snapshot number:
SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=78
USE_Z=-1  # if -1, the redshift is set by SNAPNUM

# Directory with the simulation output:
BASEDIR=/path/to/${SIMULATION}/output

# Calculate magnitudes with or without Charlot & Fall (2000) dust model:
USE_CF00=0  # 0 = no, 1 = yes

python ${CODEDIR}/galaxev_pipeline/stellar_photometrics.py \
    ${SUITE} ${BASEDIR} ${WRITEDIR} ${BC03_MODEL_DIR} \
    ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_SET}
