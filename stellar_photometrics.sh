#!/bin/bash

SUITE=IllustrisTNG
WRITEDIR=/n/hernquistfs3/vrodrigu/SyntheticImages/galaxev/IllustrisTNG/stellar_photometrics
BC03_MODEL_DIR=/n/home10/vrodrigu/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=/n/home10/vrodrigu/SyntheticImages/broadband_filters
FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt
CODEDIR=/n/home10/vrodrigu/SyntheticImages/galaxev_pipeline

for SNAPNUM in $(seq 95 99); do
  python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${SNAPNUM} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR}
done
