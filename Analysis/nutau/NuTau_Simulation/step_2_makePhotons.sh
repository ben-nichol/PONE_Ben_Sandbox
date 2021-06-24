#!/bin/bash
date
startsecond=$(date +%s)

echo "Starting the clsim job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "cvfms ran"

i3env=/home/users/akatil/software/V06-01-02/build/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/sim/clsim/step_2_clsim_setCrossE.py
echo "Will use script: " $script

MEDIUM=$1
RUNNUM=$2
LOGMINENERGY=5.0
LOGMAXENERGY=6.7

echo "Medium used will be " $MEDIUM
MEDIUMMODEL=/home/users/akatil/P-ONE/medium/STRAW_Andy_20200328_MattewEta

GCDFILE=/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz
INNAME=tauDecay_${LOGMINENERGY}_${LOGMAXENERGY}_PONE_Phase1_${RUNNUM}.i3.gz
INDIR=/data/p-one/akatil/step_1_medium_water
OUTNAME=step_2_${RUNNUM}_medium_water.i3.gz
OUTDIR=/data/p-one/akatil/step_2_medium_water

echo "INPUT FILE NAME : "$NUMEVENTS
echo "INPUT FILE DIR : "$LOGMINENERGY":"$LOGMAXENERGY
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "GCD FILE : "$GCDFILE

$i3env python $script -o ${OUTDIR}/${OUTNAME} -i ${INDIR}/${INNAME} -r ${RUNNUM} -g ${GCDFILE} -m ${MEDIUMMODEL}

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
