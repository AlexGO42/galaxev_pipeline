#!/bin/bash

#~ # Generic multi-epoch sample (rest-frame photometry)
#~ MOCK_TYPE=generic
#~ USE_Z=0.0
#~ FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt

#~ # Mock POGS sample (snapshot 131/95 from Illustris/IllustrisTNG)
#~ MOCK_TYPE=pogs
#~ USE_Z=0.0485236299818
#~ FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt

#~ # SDSS
#~ MOCK_TYPE=sdss
#~ USE_Z=0.0485236299818
#~ FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/sloan.txt

# KiDS
MOCK_TYPE=kids
USE_Z=0.1527487689024
FILENAME_FILTERS=$HOME/SyntheticImages/src/broadband_filters/vst.txt

BC03_MODEL_DIR=$HOME/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=$HOME/SyntheticImages/src/broadband_filters
CODEDIR=$HOME/SyntheticImages/src/galaxev_pipeline

for SUITE in Illustris IllustrisTNG; do
  #~ WRITEDIR=/n/hernquistfs3/vrodrigu/SyntheticImages/${MOCK_TYPE}/${SUITE}
  WRITEDIR=/u/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
  for USE_CF00 in 0 1; do
    python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR} ${USE_CF00} ${USE_Z}
  done
done
