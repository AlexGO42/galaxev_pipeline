#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J galex_xy_fof
#SBATCH -o galex_xy_fof.out
#SBATCH -e galex_xy_fof.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -t 2-00:00:00
#SBATCH -p p.48h
#SBATCH --exclusive
#SBATCH --mem=120000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# IllustrisTNG
SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=99
#USE_Z=-1  # usually set by SNAPNUM
USE_Z=0.0485236299818
NUM_RHALFS=7.5
NPIXELS=-1  # only set if NUM_RHALFS=-1

#~ # Illustris
#~ SUITE=Illustris
#~ SIMULATION=L75n1820FP
#~ SNAPNUM=131
#~ USE_Z=-1  # usually set by SNAPNUM
#~ NUM_RHALFS=5.0
#~ NPIXELS=-1  # only set if NUM_RHALFS=-1

LOG_MSTAR_MIN=9.5
MOCK_TYPE=galex  # generic, pogs, sdss, galex, kids
PROJ_KIND=xy
NGB=32
USE_FOF=1
USE_CF00=0

BASEDIR=/isaac/ptmp/gc/apillepi/sims.mirror_TNG/${SIMULATION}/output
AMDIR=/u/vrg/AngularMomentum/output/${SUITE}/${SIMULATION}
if [ ${MOCK_TYPE} = "generic" ]; then
  FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/panstarrs.txt
elif [ ${MOCK_TYPE} = "pogs" ]; then
  FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/panstarrs.txt
elif [ ${MOCK_TYPE} = "sdss" ]; then
  FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/sloan.txt
elif [ ${MOCK_TYPE} = "kids" ]; then
  FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/vst.txt
else
  FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/${MOCK_TYPE}.txt
fi
STELLAR_PHOTOMETRICS_DIR=/u/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline
NPROCESSES=${SLURM_NTASKS_PER_NODE}

python ${CODEDIR}/create_images.py \
    ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} \
    ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} \
    ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
    ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
