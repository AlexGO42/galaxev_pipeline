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

# The following parameter specifies the projection. The possible values are
# 'xy', 'yz', 'zx', 'faceon', 'edgeon', and 'planar':
PROJ_KIND=xy

# Number of neighbors that determines the adaptive smoothing scale:
NGB=32

# The field of view can be given as the number of stellar half-mass
# radii (rhalf) in each direction from the center of the image (e.g.
# NUM_RHALFS=7.5 for a total of 15 rhalfs across), as a multiple of the
# halo radius R200 (measured from the halo center), or as the *total*
# number of pixels on each side of the image. The user should only specify
# *one* of the following parameters (and set the others to -1):
NUM_RHALFS=-1
NUM_R200=-1
NPIXELS=224

# In cosmological simulations, one is often interested in generating
# images for all galaxies above a given stellar or halo mass, which can be
# specified with the parameters LOG_MSTAR_MIN and/or LOG_M200_MIN.
# Alternatively, if LOG_MSTAR_MIN=-1 and LOG_M200_MIN=-1, a custom set of
# Subfind IDs can be given in a custom text file specified by the user.
LOG_MSTAR_MIN=12.0
LOG_M200_MIN=-1
FILENAME_IDS_CUSTOM=subfind_ids.txt  # only used if the previous two are -1

# If the following is True (i.e., 1), only create images for central galaxies
# (although satellite galaxies from the same halo that happen to be in the
# field of view are also shown if the USE_FOF option is active):
CENTRALS_ONLY=0  # 0 = no, 1 = yes

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
    ${ARCSEC_PER_PIXEL} ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NUM_R200} \
    ${NPIXELS} ${LOG_MSTAR_MIN} ${LOG_M200_MIN} ${FILENAME_IDS_CUSTOM} \
    ${CENTRALS_ONLY} ${USE_FOF} ${USE_CF00} ${NPROCESSES} ${VERBOSE}
