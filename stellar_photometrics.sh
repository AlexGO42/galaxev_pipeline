#!/bin/bash

#~ # Mock POGS sample (snapshot 131/95 from Illustris/IllustrisTNG)
#~ MOCK_TYPE=pogs
#~ USE_Z=0.0485236299818

# Generic multi-epoch sample (rest-frame photometry)
MOCK_TYPE=generic
USE_Z=0.0

WRITEDIR=/n/hernquistfs3/vrodrigu/SyntheticImages/${MOCK_TYPE}
BC03_MODEL_DIR=/n/home10/vrodrigu/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=/n/home10/vrodrigu/SyntheticImages/broadband_filters
FILENAME_FILTERS=/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt
CODEDIR=/n/home10/vrodrigu/SyntheticImages/galaxev_pipeline
USE_CF00=0

python ${CODEDIR}/stellar_photometrics.py ${USE_Z} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR} ${USE_CF00}
