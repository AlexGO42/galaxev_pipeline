#!/bin/bash

# An example batch script that calls create_images.py with an
# appropriate set of parameters.

# The user-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# Root folder for this set of synthetic images:
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Specify simulation and snapshot number:
SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=78
USE_Z=-1  # if -1, the redshift is set by SNAPNUM

# Directory with the simulation output:
BASEDIR=/path/to/simulation/output

# Directory with the angular momentum vectors. This is only needed for
# face-on/edge-on projections and can be given a dummy (but not null)
# value if not needed, as we do in this example:
AMDIR=${BASEDIR}  # dummy value

# Pixel scale of the instrument or observational data of interest:
ARCSEC_PER_PIXEL=0.168  # https://hsc.mtk.nao.ac.jp/ssp/survey

# The field of view can be given either as the number of stellar half-mass
# radii (rhalf) in each direction from the center of the image (e.g.
# NUM_RHALFS=7.5 for a total of 15 rhalfs across) or as the *total*
# number of pixels on each side of the image. The user should only specify
# one of these parameters (and set the other one to -1):
NUM_RHALFS=-1  # only set if NPIXELS=-1
NPIXELS=224  # only set if NUM_RHALFS=-1

# In cosmological simulations, one is often interested in generating
# images for all galaxies above a given stellar mass, which can be
# specified with the parameter LOG_MSTAR_MIN. Alternatively, if
# LOG_MSTAR_MIN=-1, a custom set of Subfind IDs can be specified in the
# text file subfind_ids.txt, although the file has to be placed in the folder
# where the synthetic images will be created, which in this example is
# ${WRITEDIR}/IllustrisTNG/L35n2160TNG/snapnum_096/galaxev/xy/subfind_ids.txt
LOG_MSTAR_MIN=12.0  # if -1, read Subfind IDs from file

# The following parameter specifies the projection. The possible values are
# 'xy', 'yz', 'zx', 'faceon', 'edgeon', and 'planar':
PROJ_KIND=xy

# Number of neighbors that determines the adaptive smoothing scale:
NGB=32

# Whether to also show other galaxies (plus diffuse light) from the same
# parent halo (friends-of-friends group):
USE_FOF=1  # 0 = no, 1 = yes

# Use the precalculated magnitudes with or without Charlot & Fall (2000) dust model:
USE_CF00=0  # 0 = no, 1 = yes

# Number of parallel MPI processes (1 for serial execution):
NPROCESSES=1

# If verbose, print some additional output:
VERBOSE=1  # 0 = no, 1 = yes

# Create images, with or without MPI
if [ ${NPROCESSES} -eq 1 ]; then
  EXEC=python
else
  EXEC="mpiexec -n ${NPROCESSES} python"
fi
${EXEC} ${CODEDIR}/galaxev_pipeline/create_images.py \
    ${SUITE} ${SIMULATION} ${BASEDIR} ${AMDIR} \
    ${WRITEDIR} ${CODEDIR}/galaxev_pipeline ${SNAPNUM} ${USE_Z} \
    ${ARCSEC_PER_PIXEL} ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} \
    ${LOG_MSTAR_MIN} ${USE_FOF} ${USE_CF00} ${NPROCESSES} ${VERBOSE}
