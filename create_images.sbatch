#!/bin/bash

#SBATCH --mail-user=vrg@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -J pogs_TNG_cf00
#SBATCH -o pogs_TNG_cf00.out
#SBATCH -e pogs_TNG_cf00.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
#SBATCH --exclusive
#SBATCH --mem=240000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# IllustrisTNG
FILESYSTEM=hernquistfs3
SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=95

#~ # Illustris
#~ FILESYSTEM=ghernquist
#~ SUITE=Illustris
#~ SIMULATION=L75n1820FP
#~ SNAPNUM=131

MOCK_TYPE=pogs
PROJ_KIND=xy
NGB=32
NUM_RHALFS=15
USE_CF00=1

BASEDIR=/n/${FILESYSTEM}/${SUITE}/Runs/${SIMULATION}/output
AMDIR=/n/hernquistfs3/vrodrigu/AngularMomentum/output/${SUITE}/${SIMULATION}
FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt
STELLAR_PHOTOMETRICS_DIR=/n/hernquistfs3/vrodrigu/SyntheticImages/${MOCK_TYPE}/${SUITE}
WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
CODEDIR=/n/home10/vrodrigu/SyntheticImages/galaxev_pipeline
# Must match "--ntasks-per-node":
NPROCESSES=64

python ${CODEDIR}/create_images.py ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
