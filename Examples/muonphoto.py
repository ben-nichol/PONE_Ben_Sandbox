#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT combo/V00-00-04

import argparse
from os.path import expandvars
import os, sys, random
from I3Tray import *
import random
from icecube import icetray, dataclasses, dataio, simclasses
from icecube import phys_services, sim_services
from icecube import clsim

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument("-i", "--infile",default="./test_input.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-o", "--outfile",default="./test_output.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type=int, default=1, help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",type=int, default=1, help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Read in GCD file")
parser.add_argument("-e", "--efficiency", type=float,default=1.0,help="DOM Efficiency ... the same as UnshadowedFraction")
parser.add_argument("-m", "--icemodel", default="spice_3.2.1",help="Ice model (spice_mie, spice_lea, etc)")
parser.add_argument("-c", "--crossenergy", type=float,default=200.0,help="The cross energy where the hybrid clsim approach will be used")
parser.add_argument("-f", "--frames", type=int,default=100,help="N Frames")
args = parser.parse_args()
count = 0
CPU=False

photon_series = "I3Photons"
#print 'CUDA devices: ', options.DEVICE
tray = I3Tray()

# Now fire up the random number generator with that seed
#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                seed = int(args.runnumber),
               nstreams = int(4e7),
                streamnum = int(args.runnumber))

tray.context['I3RandomService'] = randomService

def PrintMessage(frame,message="") :
  print(message)
  count += 1
  return True

outfile = args.outfile + "_PhotonProp_"+str(args.runnumber)+".i3.gz"

infile = args.infile + str(args.runnumber)+".i3.gz"

icemodel_path =  args.icemodel

gcd_file = dataio.I3File(args.gcdfile)

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile,infile]
            )

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                #UseCPUs=True,
                UseGPUs=True,
                #UseOnlyDeviceNumber=[1],
                #OpenCLDeviceList=[0],
                MCTreeName="I3MCTree",
                #OutputMCTreeName="I3MCTree_clsim",
                #FlasherInfoVectName="I3FlasherInfo",
                #FlasherPulseSeriesName="PhotonBomb",
                #MMCTrackListName="MMCTrackList",
                PhotonSeriesName=photon_series,
                #ParallelEvents=1000,
                RandomService=randomService,
                IceModelLocation="/home/users/tmcelroy/pone_offline/WaterOpticalModel/STRAW_Andy_20200328_MattewEta",
                #IceModelLocation=mediumProperties,
                #UnWeightedPhotons=True, #turn off optimizations
                #UseGeant4=True,
                CrossoverEnergyEM=0.1,
                PhotonHistoryEntries=0,
                #CrossoverEnergyHadron=float(options.CROSSENERGY),
                StopDetectedPhotons=True,
                #UseHoleIceParameterization=False, # Apply it when making hits!
                #HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.flasher_p1_0.30_p2_-1"),
                DoNotParallelize=False,
                DOMOversizeFactor=1.,
                UnshadowedFraction=1., #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile=args.gcdfile #gcd_file
                )

tray.AddModule("I3Writer","writer",
#               SkipKeys=SkipKeys,
               Filename =  outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute(10)
tray.Finish()
