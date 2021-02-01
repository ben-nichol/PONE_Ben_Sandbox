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
parser.add_argument("-i", "--infile",default="./test_input.i3", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type="string", default="1", dest="RUNNUMBER", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",type="string",default="1", help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Read in GCD file")
parser.add_argument("-e","--efficiency", type="float",default=1.0,help="DOM Efficiency ... the same as UnshadowedFraction")
parser.add_argument("-m","--icemodel", default="spice_3.2.1",help="Ice model (spice_mie, spice_lea, etc)")
parser.add_argument("-c","--crossenergy", type="float",default=200.0,help="The cross energy where the hybrid clsim approach will be used")

args = parser.parse_args()

CPU=False

photon_series = "I3Photons"
#print 'CUDA devices: ', options.DEVICE
tray = I3Tray()

# Now fire up the random number generator with that seed
#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                seed = 1234*args.runnumber,
                nstreams = 1000000,
                streamnum = args.runnumber)

tray.context['I3RandomService'] = randomService

def BasicHitFilter(frame):
    hits = 0
    if frame.Has(photon_series):
       hits = len(frame.Get(photon_series))
    if hits>0:
       return True
    else:
       return False


outfile = args.outfile
infile = args.infile

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile, infile]
            )

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

icemodel_path =  args.icemodel

gcd_file = dataio.I3File(args.gcdfile)

tray.AddModule("I3InfiniteSource","streams",
               Prefix=args.gcdfile, 
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               RunNumber=args.runnumber,
               EventID=1,
               IncrementEventID=True)

tray.AddModule(PhotonBomb, "customFlasher",
               FlasherPulseSeriesName = "PhotonBomb",
               PhotonsPerPulse = 1.e3,
               )

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseCPUs=CPU,
                UseGPUs=True,
                #UseOnlyDeviceNumber=[1],
                #OpenCLDeviceList=[0],
                #MCTreeName="I3MCTree",
                OutputMCTreeName="I3MCTree_clsim",
                #FlasherInfoVectName="I3FlasherInfo",
                FlasherPulseSeriesName="PhotonBomb",
                #MMCTrackListName="MMCTrackList",
                PhotonSeriesName=photon_series,
                ParallelEvents=1000,
                RandomService=randomService,
                IceModelLocation=icemodel_path,
                #IceModelLocation=mediumProperties,
                #UnWeightedPhotons=True, #turn off optimizations
                UseGeant4=True,
                CrossoverEnergyEM=0.1,
                #PhotonHistoryEntries=1000,
                #CrossoverEnergyHadron=float(options.CROSSENERGY),
                StopDetectedPhotons=True,
                #UseHoleIceParameterization=False, # Apply it when making hits!
                #HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.flasher_p1_0.30_p2_-1"),
                DoNotParallelize=False,
                DOMOversizeFactor=1.,
                UnshadowedFraction=1., #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile=gcd_file,
                ExtraArgumentsToI3CLSimModule={
                    #"UseHardcodedDeepCoreSubdetector":True, #may save some GPU memory
                    #"EnableDoubleBuffering":True,
                    "DoublePrecision":False, #will impact performance if true
                    "StatisticsName":"clsim_stats",
                    "IgnoreDOMIDs":[],
                    #"SaveAllPhotons":True,
                    }
                )


# Tested that all frames go through CLSIM. Removing the ones without any hits to save space.
tray.AddModule(BasicHitFilter, 'FilterNullPhotons', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

SkipKeys = ["I3MCTree_bak"]

tray.AddModule("I3Writer","writer",
               SkipKeys=SkipKeys,
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
