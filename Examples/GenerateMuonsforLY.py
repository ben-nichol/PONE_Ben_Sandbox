#!/usr/bin/python

import numpy as np
import scipy.stats as st
from icecube import dataclasses, phys_services, dataio, icetray
from icecube.icetray import I3Units
import os, sys
from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                    dest="OUTFILE",help="Write output to OUTFILE (.i3{.gz} format")
parser.add_option("-f", "--frames",type=int, default=100000,dest="NUM_FRAMES",
                    help="Number of frames/MCTrees you would like to make")
parser.add_option("-r", "--radius", type=float, default="50.0", dest="RAD",
                    help="Radius to generate electrons in")
parser.add_option("-t", "--timescale",type="float",default="1000.0",dest="TIME",
                    help="Time length of 1 frame in microseconds, default 1000us")

(options,args) = parser.parse_args()
rng = np.random

num_frames = options.NUM_FRAMES 
# start frame loop, look up how to make new frames

#open new .i3 file (outfile) for writing
outfile = dataio.I3File(options.OUTFILE,'w')

for i in range(num_frames):
    frame = icetray.I3Frame('Q')

    tree = dataclasses.I3MCTree()

    theta = np.pi*2.0*rng.random_sample()
    h = 500*(-1.0 + 2.0*rng.random_sample())
    pos = [250.0*np.cos(theta),250.0*np.sin(theta),h]
    l = np.sqrt(pos[0]**2.0+pos[1]**2.0+pos[2]**2.0)
    

    muonpos = dataclasses.I3Position(pos[0],pos[1],pos[2])
    muondir = dataclasses.I3Direction(-pos[0]/l,-pos[1]/l,-pos[2]/l)

    particle = dataclasses.I3Particle()
    particle.type = dataclasses.I3Particle.MuMinus
    particle.energy = 100.0*I3Units.GeV
    particle.pos = muonpos
    particle.dir = muondir
    particle.time = 100.0*I3Units.ns
    particle.length = float('nan')
    particle.location_type_string="InIce"
    particle.shape = dataclasses.I3Particle.ParticleShape.MCTrack
    tree.add_primary(particle)
    frame.Put('I3MCTree_preMuonProp',tree)
    outfile.push(frame)

outfile.close()
