#!/bin/bash

date
startsecond=$(date +%s)

echo "Starting the NuGen job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "cvfms ran"

i3env=/home/users/akatil/software/V06-01-02/build/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/Tau_Propagation/NuTau.py
echo "Will use script: " $script

NUMEVENTS=$1
LOGMINENERGY=$2
LOGMAXENERGY=$3
RUNNUM=$4

echo "Number of events: " $NUMEVENTS

OUTNAME=tauDecay_${LOGMINENERGY}_${LOGMAXENERGY}_PONE_Phase1_${RUNNUM}
OUTDIR=/data/p-one/akatil/step_1_medium_water

echo "NUMBER OF EVENTS : "$NUMEVENTS
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "LOG ENERGY RANGE : "$LOGMINENERGY":"$LOGMAXENERGY

$i3env python $script -emin $LOGMINENERGY -emax $LOGMAXENERGY -n $NUMEVENTS -o ${OUTDIR}/${OUTNAME} -r $RUNNUM

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
