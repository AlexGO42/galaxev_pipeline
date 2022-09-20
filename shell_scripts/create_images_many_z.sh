#!/bin/bash

#SBATCH --mail-user=vrodgom.astro@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J hsc_xy
#SBATCH -o hsc_xy.out
#SBATCH -e hsc_xy.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
#SBATCH --exclusive
#SBATCH --mem=240000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=${SLURM_NTASKS_PER_NODE}

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

# Create associative array for NPIXELS such that the physical scale is > 160 kpc and NPIXELS is multiple of 32
# (https://stackoverflow.com/questions/14370133/is-there-a-way-to-create-key-value-pairs-in-bash-script/23697848):
declare -A npixels_arr
npixels_arr[0]=384
npixels_arr[1]=320
npixels_arr[2]=256  # z = 12
npixels_arr[3]=256  # z = 11
npixels_arr[4]=224  # z = 10
npixels_arr[5]=224
npixels_arr[6]=224  # z = 9
npixels_arr[7]=224
npixels_arr[8]=224  # z = 8
npixels_arr[9]=192
npixels_arr[10]=192
npixels_arr[11]=192  # z = 7
npixels_arr[12]=192
npixels_arr[13]=192  # z = 6
npixels_arr[14]=192
npixels_arr[15]=160
npixels_arr[16]=160
npixels_arr[17]=160  # z = 5
npixels_arr[18]=160
npixels_arr[19]=160
npixels_arr[20]=160
npixels_arr[21]=160  # z = 4
npixels_arr[22]=160
npixels_arr[23]=128
npixels_arr[24]=128
npixels_arr[25]=128  # z = 3
npixels_arr[26]=128
npixels_arr[27]=128
npixels_arr[28]=128
npixels_arr[29]=128
npixels_arr[30]=128
npixels_arr[31]=128
npixels_arr[32]=128
npixels_arr[33]=128  # z = 2
npixels_arr[34]=128
npixels_arr[35]=128
npixels_arr[36]=128
npixels_arr[37]=128
npixels_arr[38]=128
npixels_arr[39]=128
npixels_arr[40]=128  # z = 1.5
npixels_arr[41]=128
npixels_arr[42]=128
npixels_arr[43]=128
npixels_arr[44]=128
npixels_arr[45]=128
npixels_arr[46]=128
npixels_arr[47]=128
npixels_arr[48]=128
npixels_arr[49]=128
npixels_arr[50]=128  # z = 1.0
npixels_arr[51]=128
npixels_arr[52]=128
npixels_arr[53]=128
npixels_arr[54]=128
npixels_arr[55]=128
npixels_arr[56]=128
npixels_arr[57]=128
npixels_arr[58]=128
npixels_arr[59]=160  # z = 0.7
npixels_arr[60]=160
npixels_arr[61]=160
npixels_arr[62]=160
npixels_arr[63]=160
npixels_arr[64]=160
npixels_arr[65]=160
npixels_arr[66]=160
npixels_arr[67]=160  # z = 0.5
npixels_arr[68]=160
npixels_arr[69]=160
npixels_arr[70]=192
npixels_arr[71]=192
npixels_arr[72]=192  # z = 0.4
npixels_arr[73]=192
npixels_arr[74]=192
npixels_arr[75]=192
npixels_arr[76]=224
npixels_arr[77]=224
npixels_arr[78]=224  # z = 0.3
npixels_arr[79]=224
npixels_arr[80]=256
npixels_arr[81]=256
npixels_arr[82]=256
npixels_arr[83]=288
npixels_arr[84]=288  # z = 0.2
npixels_arr[85]=320
npixels_arr[86]=320
npixels_arr[87]=352
npixels_arr[88]=384
npixels_arr[89]=416
npixels_arr[90]=480
npixels_arr[91]=512  # z = 0.1
npixels_arr[92]=608
npixels_arr[93]=672
npixels_arr[94]=832
npixels_arr[95]=992
npixels_arr[96]=1376
npixels_arr[97]=1920
npixels_arr[98]=4736
npixels_arr[99]=199903361914974976  # z = 0.0 (duh)

BASEDIR=/n/holylfs05/LABS/hernquist_lab/${SUITE}/Runs/${SIMULATION}/output
AMDIR=/n/holystore01/LABS/hernquist_lab/Lab/vrodrigu/AngularMomentum/output/${SUITE}/${SIMULATION}
FILTER_DIR=${HOME}/SyntheticImages/src/broadband_filters
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline
NPROCESSES=${SLURM_NTASKS_PER_NODE}

#for MOCK_TYPE in candels_wfc3 candels_acs; do
for MOCK_TYPE in hsc; do
  #FILENAME_FILTERS=${FILTER_DIR}/${MOCK_TYPE}.txt
  FILENAME_FILTERS=${FILTER_DIR}/subaru.txt
  STELLAR_PHOTOMETRICS_DIR=/n/holyscratch01/hernquist_lab/Lab/vrodrigu/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
  WRITEDIR=${STELLAR_PHOTOMETRICS_DIR}/${SIMULATION}
  for SNAPNUM in $(seq 40 91); do
    NPIXELS=${npixels_arr[${SNAPNUM}]}
    mpiexec -n ${NPROCESSES} python ${CODEDIR}/create_images.py \
        ${SUITE} ${BASEDIR} ${AMDIR} ${FILENAME_FILTERS} \
        ${STELLAR_PHOTOMETRICS_DIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} ${USE_Z} \
        ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${NPIXELS} ${LOG_MSTAR_MIN} ${USE_FOF} \
        ${USE_CF00} ${MOCK_TYPE} ${NPROCESSES}
    echo "Finished for snapshot ${SNAPNUM}."
  done
done
