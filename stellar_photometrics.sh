#!/bin/bash

# Create any name for a new set of synthetic images:
MOCK_SET=hsc

# A folder corresponding to the above set of synthetic images
# (can include many simulations, snapshots, etc.)
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Simulation snapshot, although not really needed yet (only for the filename):
SNAPNUM=91

# Observation redshift, which usually corresponds to snapnum:
USE_Z=0.0994018026302

# Directory with the GALAXEV model files:
BC03_MODEL_DIR=/path/to/bc03/Padova1994/chabrier

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Illustris or IllustrisTNG (other simulations eventually):
SUITE=IllustrisTNG

# Calculate magnitudes with and without Charlot & Fall (2000) dust model.
for USE_CF00 in 0 1; do
  python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} \
         ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_SET}
done
