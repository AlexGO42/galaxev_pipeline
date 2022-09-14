#!/bin/bash

# This example is similar to stellar_photometrics.sbatch but calculates
# apparent magnitudes for a large number of redshifts, corresponding to
# all IllustrisTNG snapshots between 40 and 91 (z = 1.5 to z = 0.1).

# Associative array (only works with Bash > 4.0) with all the
# IllustrisTNG snapshot redshifts:
declare -A redshifts
redshifts[0]=2.00464909888e+01
redshifts[1]=1.49891732400e+01
redshifts[2]=1.19802133153e+01
redshifts[3]=1.09756432941e+01
redshifts[4]=9.99659046619e+00
redshifts[5]=9.38877127194e+00
redshifts[6]=9.00233985416e+00
redshifts[7]=8.44947629437e+00
redshifts[8]=8.01217294887e+00
redshifts[9]=7.59510714987e+00
redshifts[10]=7.23627606617e+00
redshifts[11]=7.00541704554e+00
redshifts[12]=6.49159774567e+00
redshifts[13]=6.01075739884e+00
redshifts[14]=5.84661374788e+00
redshifts[15]=5.52976580795e+00
redshifts[16]=5.22758097313e+00
redshifts[17]=4.99593346816e+00
redshifts[18]=4.66451770247e+00
redshifts[19]=4.42803373661e+00
redshifts[20]=4.17683491473e+00
redshifts[21]=4.00794511147e+00
redshifts[22]=3.70877426464e+00
redshifts[23]=3.49086136926e+00
redshifts[24]=3.28303305796e+00
redshifts[25]=3.00813107163e+00
redshifts[26]=2.89578500573e+00
redshifts[27]=2.73314261732e+00
redshifts[28]=2.57729027160e+00
redshifts[29]=2.44422570455e+00
redshifts[30]=2.31611074396e+00
redshifts[31]=2.20792547238e+00
redshifts[32]=2.10326965260e+00
redshifts[33]=2.00202813925e+00
redshifts[34]=1.90408954353e+00
redshifts[35]=1.82268925262e+00
redshifts[36]=1.74357057433e+00
redshifts[37]=1.66666955611e+00
redshifts[38]=1.60423452207e+00
redshifts[39]=1.53123902916e+00
redshifts[40]=1.49551216650e+00
redshifts[41]=1.41409822037e+00
redshifts[42]=1.35757666740e+00
redshifts[43]=1.30237845991e+00
redshifts[44]=1.24847261425e+00
redshifts[45]=1.20625808078e+00
redshifts[46]=1.15460271236e+00
redshifts[47]=1.11415056377e+00
redshifts[48]=1.07445789455e+00
redshifts[49]=1.03551044566e+00
redshifts[50]=9.97294225782e-01
redshifts[51]=9.50531351585e-01
redshifts[52]=9.23000816178e-01
redshifts[53]=8.86896937575e-01
redshifts[54]=8.51470900625e-01
redshifts[55]=8.16709979012e-01
redshifts[56]=7.91068248946e-01
redshifts[57]=7.57441372616e-01
redshifts[58]=7.32636182022e-01
redshifts[59]=7.00106353719e-01
redshifts[60]=6.76110411213e-01
redshifts[61]=6.44641840685e-01
redshifts[62]=6.21428745243e-01
redshifts[63]=5.98543288188e-01
redshifts[64]=5.75980845108e-01
redshifts[65]=5.46392183141e-01
redshifts[66]=5.24565820434e-01
redshifts[67]=5.03047523245e-01
redshifts[68]=4.81832943421e-01
redshifts[69]=4.60917794181e-01
redshifts[70]=4.40297849248e-01
redshifts[71]=4.19968941997e-01
redshifts[72]=3.99926964614e-01
redshifts[73]=3.80167867260e-01
redshifts[74]=3.60687657262e-01
redshifts[75]=3.47853841858e-01
redshifts[76]=3.28829724206e-01
redshifts[77]=3.10074120128e-01
redshifts[78]=2.97717684517e-01
redshifts[79]=2.73353346578e-01
redshifts[80]=2.61343256161e-01
redshifts[81]=2.43540181555e-01
redshifts[82]=2.25988386260e-01
redshifts[83]=2.14425035514e-01
redshifts[84]=1.97284182376e-01
redshifts[85]=1.80385261706e-01
redshifts[86]=1.69252033244e-01
redshifts[87]=1.52748768902e-01
redshifts[88]=1.41876203970e-01
redshifts[89]=1.25759332411e-01
redshifts[90]=1.09869940459e-01
redshifts[91]=9.94018026302e-02
redshifts[92]=8.38844307975e-02
redshifts[93]=7.36613846564e-02
redshifts[94]=5.85073227945e-02
redshifts[95]=4.85236299818e-02
redshifts[96]=3.37243718735e-02
redshifts[97]=2.39744283828e-02
redshifts[98]=9.52166696794e-03
redshifts[99]=2.22044604925e-16

# Any user-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# A folder to store all the output corresponding to this set of synthetic
# images (with subdirectories corresponding to different simulations,
# snapshots, etc., as will become clear later):
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Directory with the GALAXEV model files:
BC03_MODEL_DIR=/path/to/bc03/Padova1994/chabrier

# Directory with the galaxev_pipeline code:
CODEDIR=/path/to/galaxev_pipeline

# Illustris or IllustrisTNG (other simulations eventually):
SUITE=IllustrisTNG

# Calculate magnitudes with or without Charlot & Fall (2000) dust model:
USE_CF00=0  # 0 = no, 1 = yes

for SNAPNUM in $(seq 40 91); do
  USE_Z=${redshifts[${SNAPNUM}]}
  python ${CODEDIR}/stellar_photometrics.py ${SUITE} ${WRITEDIR} ${BC03_MODEL_DIR} \
       ${USE_CF00} ${SNAPNUM} ${USE_Z} ${MOCK_SET}
  echo "Finished for snapshot ${SNAPNUM}."
done
