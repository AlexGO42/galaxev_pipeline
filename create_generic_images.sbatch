#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J generic_xy
#SBATCH -o generic_xy.out
#SBATCH -e generic_xy.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -t 2-00:00:00
#SBATCH -p p.48h
#SBATCH --exclusive
#SBATCH --mem=120000

# OMP_NUM_THREADS controls the number of threads your application uses:
#export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
OMP_NUM_THREADS=16

# IllustrisTNG
SUITE=IllustrisTNG
SIMULATION=L75n1820TNG
SNAPNUM=99
NUM_RHALFS=7.5
NPIXELS=-1  # only set if NUM_RHALFS=-1

LOG_MSTAR_MIN=9.0  # if -1, read Subfind IDs from file
PROJ_KIND=xy
NGB=32
USE_FOF=0

BASEDIR=/isaac/ptmp/gc/apillepi/sims.mirror_TNG/${SIMULATION}/output
AMDIR=/u/vrg/AngularMomentum/output/${SUITE}/${SIMULATION}
WRITEDIR=/isaac/ptmp/gc/vrg/SyntheticImages/output/generic/${SUITE}/${SIMULATION}
CODEDIR=${HOME}/SyntheticImages/src/galaxev_pipeline
#NPROCESSES=${SLURM_NTASKS_PER_NODE}
NPROCESSES=16

# mpiexec does not seem to obey the stop signal
# (https://bugzilla.redhat.com/show_bug.cgi?id=1151657).
# The following should give me time to send the job to the background
# before mpiexec is called:
echo "Sleeping for 10 seconds..."
sleep 10

mpiexec -n ${NPROCESSES} python ${CODEDIR}/create_generic_images.py \
    ${SUITE} ${SIMULATION} ${BASEDIR} ${AMDIR} ${WRITEDIR} ${CODEDIR} ${SNAPNUM} \
    ${PROJ_KIND} ${NGB} ${NUM_RHALFS} ${LOG_MSTAR_MIN} ${USE_FOF} \
    ${NPROCESSES}
