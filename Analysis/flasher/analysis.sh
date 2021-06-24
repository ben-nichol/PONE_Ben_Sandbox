#!/bin/bash

date
startsecond=$(date +%s)

echo "Starting the analysis job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "cvfms ran"

i3env=/home/users/akatil/software/V06-01-02/build/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/medium/runCustomFlashes.py
echo "Will use script: " $script

NUMFRAMES=$1
RUNNUM=$2

BLUE=1

script2=/home/users/akatil/P-ONE/medium/makePhotonsFromFlashes.py
echo "Will use script: " $script2

MEDIUMMODEL=/home/users/akatil/P-ONE/medium/STRAW_Andy_20200328_MattewEta
GCDFILE=/home/users/akatil/P-ONE/GCD_files/STRAW_DOM_at_zero.i3.gz

OUTNAMEF53=genUvIsotropicFlashes_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_53.i3.gz
OUTDIRF53=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta
OUTNAMEP53=genUvPhotons_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_53.i3.gz
OUTDIRP53=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta

echo "OUTPUT FILE NAME : "$OUTNAMEP53
echo "OUTPUT FILE DIR  : "$OUTDIRP53
echo "GCD FILE : "$GCDFILE

$i3env python $script2 -o ${OUTDIRP53}/${OUTNAMEP53} -i ${OUTDIRF53}/${OUTNAMEF53} -r ${RUNNUM} -g ${GCDFILE} -m ${MEDIUMMODEL}

OUTNAMEF70=genUvIsotropicFlashes_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_70.i3.gz
OUTDIRF70=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta
OUTNAMEP70=genUvPhotons_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_70.i3.gz
OUTDIRP70=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta

echo "OUTPUT FILE NAME : "$OUTNAMEP70
echo "OUTPUT FILE DIR  : "$OUTDIRP70
echo "GCD FILE : "$GCDFILE

$i3env python $script2 -o ${OUTDIRP70}/${OUTNAMEP70} -i ${OUTDIRF70}/${OUTNAMEF70} -r ${RUNNUM} -g ${GCDFILE} -m ${MEDIUMMODEL}

OUTNAMEF88=genUvIsotropicFlashes_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_88.i3.gz
OUTDIRF88=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta
OUTNAMEP88=genUvPhotons_${RUNNUM}_${BLUE}_500000Loop_2e6_spread_88.i3.gz
OUTDIRP88=/data/p-one/akatil/timeResiduals/Wavelength_20200326/UV_20200328_MattewEta

echo "OUTPUT FILE NAME : "$OUTNAMEP88
echo "OUTPUT FILE DIR  : "$OUTDIRP88
echo "GCD FILE : "$GCDFILE

$i3env python $script2 -o ${OUTDIRP88}/${OUTNAMEP88} -i ${OUTDIRF88}/${OUTNAMEF88} -r ${RUNNUM} -g ${GCDFILE} -m ${MEDIUMMODEL}

echo "2nd Step Done"

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
