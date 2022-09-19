#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J kids_xy_fof
#SBATCH -o kids_xy_fof.out
#SBATCH -e kids_xy_fof.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -t 2-00:00:00
#SBATCH -p p.48h
#SBATCH --exclusive
#SBATCH --mem=120000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# User-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# Root folder to store the corresponding output:
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Directory with the simulation output
BASEDIR=/path/to/simulation/output

# Directory with the angular momentum vectors. This is only needed for
# face-on/edge-on projections and can be given a dummy value if not needed:
AMDIR=${BASEDIR}

SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=91
USE_Z=-1  # usually set by SNAPNUM
ARCSEC_PER_PIXEL=0.168  # https://hsc.mtk.nao.ac.jp/ssp/survey
NUM_RHALFS=-1
NPIXELS=128  # only set if NUM_RHALFS=-1

LOG_MSTAR_MIN=9.0  # if -1, read Subfind IDs from file
PROJ_KIND=xy
NGB=32
USE_FOF=1
USE_CF00=0

NPROCESSES=${SLURM_NTASKS_PER_NODE}

mpiexec -n ${NPROCESSES} python ${CODEDIR}/create_images.py \
    ${SUITE} ${SIMULATION} ${BASEDIR} ${AMDIR} \
    ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} ${ARCSEC_PER_PIXEL} \
    ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
    ${USE_CF00} ${NPROCESSES}
