#!/bin/bash

# For HST (CANDELS) and HSC, we work with the following snapshots:
####################
# snapnum, redshift
# 25, 3.00813107163
# 29, 2.44422570455
# 33, 2.00202813925
# 40, 1.4955121665
# 50, 0.997294225782
# 67, 0.503047523245
####################

# Create associative array (only works with Bash > 4.0)
# (https://stackoverflow.com/questions/14370133/is-there-a-way-to-create-key-value-pairs-in-bash-script/23697848):
declare -A redshifts
redshifts[25]=3.00813107163
redshifts[29]=2.44422570455
redshifts[33]=2.00202813925
redshifts[40]=1.4955121665
redshifts[50]=0.997294225782
redshifts[67]=0.503047523245

BC03_MODEL_DIR=${HOME}/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=${HOME}/SyntheticImages/src/broadband_filters
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline

SUITE=IllustrisTNG
USE_CF00=0
#for MOCK_TYPE in candels_wfc3 candels_acs; do
for MOCK_TYPE in hsc; do
  WRITEDIR=/isaac/ptmp/gc/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
  #FILENAME_FILTERS=${FILTER_DIR}/${MOCK_TYPE}.txt
  FILENAME_FILTERS=${FILTER_DIR}/subaru.txt
  for SNAPNUM in 25 29 33 40 50 67; do
    USE_Z=${redshifts[${SNAPNUM}]}
    python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR} ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_TYPE}
    echo "Finished for snapshot ${SNAPNUM}."
  done
done
