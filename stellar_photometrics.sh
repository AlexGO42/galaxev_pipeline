#!/bin/bash

## Generic multi-epoch sample (rest-frame photometry)
#MOCK_TYPE=generic
#SNAPNUM=99  # just for the filename
#USE_Z=0.0
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/panstarrs.txt

## Mock POGS sample (snapshot 131/95 from Illustris/IllustrisTNG)
#MOCK_TYPE=pogs
#SNAPNUM=95  # just for the filename
#USE_Z=0.0485236299818  # must be set by hand
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/panstarrs.txt

## SDSS
#MOCK_TYPE=sdss
#SNAPNUM=95  # just for the filename
#USE_Z=0.0485236299818  # must be set by hand
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/sloan.txt

## SDSS BCG
#MOCK_TYPE=sdss_bcg
#SNAPNUM=91  # just for the filename
#USE_Z=0.0994018026302  # must be set by hand
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/sloan.txt

## GALEX
#MOCK_TYPE=galex
#SNAPNUM=95  # just for the filename
#USE_Z=0.0485236299818  # must be set by hand
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/galex.txt

## KiDS Lingyu
#MOCK_TYPE=kids
#SNAPNUM=87  # just for the filename
#USE_Z=0.1527487689024  # must be set by hand
#FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/vst.txt

# KiDS
MOCK_TYPE=kids
SNAPNUM=91  # just for the filename
USE_Z=0.0994018026302  # must be set by hand
FILENAME_FILTERS=${HOME}/SyntheticImages/src/broadband_filters/kids.txt

BC03_MODEL_DIR=${HOME}/galaxev_code/bc03/Padova1994/chabrier
FILTER_DIR=${HOME}/SyntheticImages/src/broadband_filters
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline

SUITE=IllustrisTNG
WRITEDIR=/isaac/ptmp/gc/vrg/SyntheticImages/output/${MOCK_TYPE}/${SUITE}
for USE_CF00 in 0 1; do
  python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} ${FILTER_DIR} ${FILENAME_FILTERS} ${CODEDIR} ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_TYPE}
done
