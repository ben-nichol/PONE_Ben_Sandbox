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
script=/home/users/akatil/P-ONE/sim/clsim/step3_mcpeHits.py
echo "Will use script: " $script

MEDIUM=$1
RUNNUM=$2

echo "Medium used will be " $MEDIUM
MEDIUMMODEL=/home/users/akatil/P-ONE/medium/STRAW_Andy_20200328_MattewEta

GCDFILE=/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz
INNAME=step_2_${RUNNUM}_medium_water.i3.gz
INDIR=/data/p-one/akatil/step_2_medium_water
OUTNAME=step_3_${RUNNUM}_medium_water_icecube.i3.gz
OUTDIR=/data/p-one/akatil/step_3_medium_water/IceCube

echo "RUNNUM : "$RUNNUM
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "GCD FILE : "$GCDFILE

$i3env python $script -o ${OUTDIR}/${OUTNAME} -i ${INDIR}/${INNAME} -r ${RUNNUM} -g ${GCDFILE} -m ${MEDIUMMODEL}

script2=/home/users/akatil/P-ONE/sim/clsim/genHitsFromI3Photons.py

GCDTYPE=PONE
DOMTYPE=IceCube
HITTHRESH=1
DOMTHRESH=1
OUTNAME=step_3_${RUNNUM}_medium_water_custom.i3.gz
OUTDIR=/data/p-one/akatil/step_3_medium_water/Custom

echo "RUNNUM : "$RUNNUM
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "GCD FILE : "$GCDFILE

$i3env python $script2 -n ${RUNNUM} -g ${GCDTYPE} -d ${DOMTYPE} -H ${HITTHRESH} -D ${DOMTHRESH} -i ${INDIR}/${INNAME} -o ${OUTDIR}/${OUTNAME}

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
