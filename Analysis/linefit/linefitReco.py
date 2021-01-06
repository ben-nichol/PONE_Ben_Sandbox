#!/usr/bin/env python


# Modified linefit file written originally by Dvir. Modified to work in my directory with some quality of life improvements. 

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
#import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
from dvir_Analysis import SimAnalysis
import argparse

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
#parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
#parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-H', '--hitThresh', default = 1, help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-D', '--domThresh', default = 3, help = "threshold of hit DOMs for the frame to be considered")
parser.add_argument( '-R', '--maxResidual', default = 100 , help = "maximum time residual allowed for the hit to be considered")
parser.add_argument( '-g', '--GCDType', help = "type of geometry used for the simulation set")
#parser.add_argument( '-i', '--input', help = "i3 file to be read in")
args = parser.parse_args()
'''
if args.GCDType == 'testString':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/testStrings/HorizTestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = '/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
elif args.GCDType == 'cube':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")
'''
gcdPath = '/home/users/ghuman/ICECUBE/gcdfiles/pentagon10040.i3'
'''
infileList = []
for i in range(int(args.minFileNum), int(args.maxFileNum)+1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/' + str(args.GCDType) + '/NuGen_step3_' + str(args.GCDType) + '_' + str(i) + '.i3.gz')
    infileList.append(infile)
'''
infileList = [dataio.I3File('/home/users/ghuman/ICECUBE/step3/NuMu/NuMu_F_pentagon10040_141_step3.zst')]
outfile = dataio.I3File('/home/users/ghuman/simAnalysis/output/I3Files/linefitReco/'+ str(args.GCDType) + '/NuGen_linefitReco_' + str(args.GCDType) + '_141.i3.gz', 'w')
gcdfile = dataio.I3File(gcdPath)
geometry = gcdfile.pop_frame()["I3Geometry"]
i = 0
for infile in infileList:
    for frame in infile:
        if SimAnalysis.passFrame(frame, geometry.omgeo.keys(), int(args.hitThresh), int(args.domThresh), float(args.maxResidual), geometry.omgeo):
            frame = SimAnalysis.writeSigHitsMapToFrame(frame, geometry.omgeo.keys(), int(args.hitThresh), int(args.domThresh), float(args.maxResidual), geometry.omgeo)
            datapoints = SimAnalysis.getLinefitDataPoints(frame, geometry)
            direction, speed, vertex = SimAnalysis.linefitParticleParams(datapoints)

            recoParticle = I3Particle()
            recoParticle.shape = I3Particle.InfiniteTrack
            recoParticle.fit_status = I3Particle.OK
            recoParticle.dir = direction
            recoParticle.speed = speed
            recoParticle.pos = vertex
            recoParticle.time = 0

            frame.Put("LineFitRecoParticle", recoParticle)
            outfile.push(frame)

outfile.close()
