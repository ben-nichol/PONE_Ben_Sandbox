#!/bin/bash

#eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
#i3env=/home/users/mens/software/V06-01-02/build/env-shell.sh
i3env=/data/p-one/jstacho/pone_offline/env-shell_Container.sh
#i3env=/home/users/tmcelroy/pone_offline/env-shell_Container.sh

script=/data/p-one/jstacho/noise/correlatedNoise/noiseGenerator/noiseModule/poneoffline/IcetrayTest.py

RUNNUM=$1
let OUTNUM=($RUNNUM)

INDIR=/data2/p-one/jstacho/noise/clsim/OutGenFast
INNAME=genHits${OUTNUM}.i3.gz
#OUTDIR=/data/p-one/jstacho/noise/correlatedNoise/noiseGenerator/noiseModule
OUTDIR=/data/p-one/jstacho/noise/clsim/genModuleOut-PONECoincs
OUTNAME=outGen${OUTNUM}.i3.gz

$i3env python $script -i ${INDIR}/${INNAME} -o ${OUTDIR}/${OUTNAME}
