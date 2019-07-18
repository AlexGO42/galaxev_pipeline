#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J sdss
#SBATCH -o sdss.out
#SBATCH -e sdss.err
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
#~ SNAPNUM=99
SNAPNUM=87
USE_Z=-1  # usually set by SNAPNUM
#~ NUM_RHALFS=7.5
NUM_RHALFS=-1
NPIXELS=128

#~ # Illustris
#~ SUITE=Illustris
#~ SIMULATION=L75n1820FP
#~ SNAPNUM=135
#~ USE_Z=-1  # usually set by SNAPNUM
#~ NUM_RHALFS=5.0
#~ NPIXELS=-1

LOG_MSTAR_MIN=9.0
MOCK_TYPE=sdss  # pogs, sdss, generic
PROJ_KIND=xy
NGB=32
USE_FOF=0
USE_CF00=0

BASEDIR=/isaac/ptmp/gc/apillepi/sims.mirror_TNG/${SIMULATION}/output
AMDIR=/u/vrg/AngularMomentum/output/${SUITE}/${SIMULATION}
if [ ${MOCK_TYPE} = "sdss" ]; then
  FILENAME_FILTERS=/u/vrg/SyntheticImages/src/broadband_filters/sloan.txt
else
  FILENAME_FILTERS=/u/vrg/SyntheticImages/src/broadband_filters/panstarrs.txt
fi
STELLAR_PHOTOMETRICS_DIR=/u/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
CODEDIR=/u/vrg/SyntheticImages/src/galaxev_pipeline
NPROCESSES=${SLURM_NTASKS_PER_NODE}

python ${CODEDIR}/create_images.py \
    ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} \
    ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} \
    ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
    ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
