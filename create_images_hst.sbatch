#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J candels_xy
#SBATCH -o candels_xy.out
#SBATCH -e candels_xy.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -t 2-00:00:00
#SBATCH -p p.48h
#SBATCH --exclusive
#SBATCH --mem=120000

# OMP_NUM_THREADS controls the number of threads your application uses:
#export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
export OMP_NUM_THREADS=16

# IllustrisTNG
SUITE=IllustrisTNG
SIMULATION=L35n2160TNG
USE_Z=-1  # set by SNAPNUM below
NUM_RHALFS=7.5
NPIXELS=-1  # only set if NUM_RHALFS=-1
LOG_MSTAR_MIN=9.0
PROJ_KIND=xy
NGB=32
USE_FOF=0
USE_CF00=0

BASEDIR=/virgo/simulations/${SUITE}/${SIMULATION}/output
AMDIR=/u/vrg/AngularMomentum/output/${SUITE}/${SIMULATION}
FILTER_DIR=${HOME}/SyntheticImages/src/broadband_filters
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline
#NPROCESSES=${SLURM_NTASKS_PER_NODE}
NPROCESSES=16

for MOCK_TYPE in candels_wfc3 candels_acs; do
  FILENAME_FILTERS=${FILTER_DIR}/${MOCK_TYPE}.txt
  STELLAR_PHOTOMETRICS_DIR=/isaac/ptmp/gc/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
  WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
  for SNAPNUM in 25 29 33 40 50 67; do
    mpiexec -n ${NPROCESSES} python ${CODEDIR}/create_images.py \
        ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} \
        ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} \
        ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
        ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
    echo "Finished for snapshot ${SNAPNUM}."
  done
done
