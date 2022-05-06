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
from Flasher.genPhotonBomb import PhotonBomb

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
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

outfile = args.outfile+"_"+str(args.runnumber)+".i3.gz"

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

def DeleteFrames(frame) :
  if frame.Has("PhotonBomb") :
    frame.Delete("PhotonBomb")
  return True

def BasicHitFilter(frame):
    hits = 0
    if frame.Has(photon_series):
       hits = len(frame[photon_series])
    if hits>0:
       return True
    else:
       return False

outfile = ""
if args.outfile[-1] == "/" :
  outfile = args.outfile + "PhotonBomb_"+str(args.runnumber)+".i3.gz"
else :
  outfile = args.outfile + "/PhotonBomb_"+str(args.runnumber)+".i3.gz" 

icemodel_path =  args.icemodel

gcd_file = dataio.I3File(args.gcdfile)

tray.Add("I3InfiniteSource", prefix = args.gcdfile)

tray.Add("I3MCEventHeaderGenerator",
             EventID=0,
             RunNumber=args.runnumber,
             IncrementEventID=True)

tray.AddModule(PhotonBomb, "customFlasher",
               FlasherPulseSeriesName = "PhotonBomb",
               PhotonsPerPulse = int(1e6),
               RandomService = randomService,
               NumPulses = int(1e3)
               )
tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseGPUs=True,
                MCTreeName="PhotonBomb",
                UseI3PropagatorService=False,
                FlasherPulseSeriesName="PhotonBomb",
                #MMCTrackListName="MMCTrackList",
                PhotonSeriesName=photon_series,
                MCPESeriesName='',
                RandomService=randomService,
                IceModelLocation=Medium.MakePoneMediumProperties(),
                UnWeightedPhotons=False,
                DOMRadius = (17.0*2.54*0.01/2.0)*icetray.I3Units.m,
                UseGeant4=False,
                CrossoverEnergyEM=None,
                PhotonHistoryEntries=0,
                StopDetectedPhotons=True,
                DoNotParallelize=False,
                WavelengthAcceptance = dom_properties.GetCLSimQETable( factor=dom_properties.GetMaxAngularAcceptance()*1.05 ),
                DOMOversizeFactor=1.0, #(17./13.),
                UnshadowedFraction=1., #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile= args.gcdfile, #gcd_file,
                )

# Tested that all frames go through CLSIM. Removing the ones without any hits to save space.
#tray.AddModule(BasicHitFilter, 'FilterNullPhotons', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

#SkipKeys = ["I3MCTree_bak"]

tray.AddModule(DeleteFrames, "DeleteStuff",Streams = [icetray.I3Frame.DAQ,
  icetray.I3Frame.Physics])

tray.AddModule("I3Writer","writer",
#               SkipKeys=SkipKeys,
               Filename =  outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute(args.frames)
tray.Finish()
