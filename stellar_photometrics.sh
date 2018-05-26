#!/bin/bash

#~ # IllustrisTNG
#~ SUITE=IllustrisTNG
#~ SNAPNUM_FIRST=95
#~ SNAPNUM_LAST=99

# Illustris
SUITE=Illustris
SNAPNUM_FIRST=131
SNAPNUM_LAST=135

WRITEDIR=/n/hernquistfs3/vrodrigu/SyntheticImages/galaxev/${SUITE}/stellar_photometrics
BC03_MODEL_DIR=/n/home10/vrodrigu/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=/n/home10/vrodrigu/SyntheticImages/broadband_filters
FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt
CODEDIR=/n/home10/vrodrigu/SyntheticImages/galaxev_pipeline
USE_CF00=1

for SNAPNUM in $(seq ${SNAPNUM_FIRST} ${SNAPNUM_LAST}); do
  python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${SNAPNUM} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR} ${USE_CF00}
done
