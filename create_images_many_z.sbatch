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
SIMULATION=L75n1820TNG
USE_Z=-1  # set by SNAPNUM below
NUM_RHALFS=-1
#~ NPIXELS=128  # only set if NUM_RHALFS=-1
LOG_MSTAR_MIN=9.0  # if -1, read Subfind IDs from file
PROJ_KIND=xy
NGB=32
USE_FOF=1
USE_CF00=0

# Create associative array for NPIXELS such that the physical scale is > 160 kpc
# (https://stackoverflow.com/questions/14370133/is-there-a-way-to-create-key-value-pairs-in-bash-script/23697848):
declare -A npixels_arr
npixels_arr[40]=128  # z = 1.5
npixels_arr[50]=128  # z = 1.0
npixels_arr[59]=160  # z = 0.7
npixels_arr[67]=160  # z = 0.5
npixels_arr[72]=192  # z = 0.4
npixels_arr[78]=256  # z = 0.3
npixels_arr[84]=320  # z = 0.2
npixels_arr[91]=512  # z = 0.1

BASEDIR=/virgotng/universe/${SUITE}/${SIMULATION}/output
AMDIR=${HOME}/AngularMomentum/output/${SUITE}/${SIMULATION}
FILTER_DIR=${HOME}/SyntheticImages/src/broadband_filters
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline
#NPROCESSES=${SLURM_NTASKS_PER_NODE}
NPROCESSES=16

#for MOCK_TYPE in candels_wfc3 candels_acs; do
for MOCK_TYPE in hsc; do
  #FILENAME_FILTERS=${FILTER_DIR}/${MOCK_TYPE}.txt
  FILENAME_FILTERS=${FILTER_DIR}/subaru.txt
  STELLAR_PHOTOMETRICS_DIR=/isaac/ptmp/gc/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
  WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
  for SNAPNUM in 40 50 59 67 72 78 84 91; do
    NPIXELS=${npixels_arr[${SNAPNUM}]}
    mpiexec -n ${NPROCESSES} python ${CODEDIR}/create_images.py \
        ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} \
        ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} \
        ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
        ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
    echo "Finished for snapshot ${SNAPNUM}."
  done
done
