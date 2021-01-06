#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars
import os, sys, random

# This script will perform a hybridCLSim propagation.
#
# NOTE: There is no bad_dom_cleaning!!!
#       This you still have to do after the propagation!!!

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="string", default="1",
                  dest="RUNNUMBER", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_option("-l", "--filenr",type="string",default="1",
                   dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_option("-g", "--gcdfile", default=os.getenv('GCDfile'),
		  dest="GCDFILE", help="Read in GCD file")
parser.add_option("-e","--efficiency", type="float",default=1.2, # Using efficiency > 1 as default so we can support systematics sets
                  dest="EFFICIENCY",help="DOM Efficiency ... the same as UnshadowedFraction")
parser.add_option("-m","--icemodel", default="spice_3.2.1",
                  dest="ICEMODEL",help="Ice model (spice_mie, spice_lea, etc)")
parser.add_option("-c","--crossenergy", type="float",default=200.0,
                  dest="CROSSENERGY",help="The cross energy where the hybrid clsim approach will be used")
parser.add_option("-t", action="store_true",  dest="GPU", default=True ,help="Run on GPUs or CPUs")


(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

options.FILENR=int(options.FILENR)
options.RUNNUMBER=int(options.RUNNUMBER)
if options.GPU:
        CPU=False
else:
        CPU=True

from I3Tray import *
import random
from icecube import icetray, dataclasses, dataio, simclasses
from icecube import phys_services, sim_services
#from icecube import diplopia
from icecube import clsim

photon_series = "I3Photons"
#print 'CUDA devices: ', options.DEVICE
tray = I3Tray()
print 'Using RUNNUMBER: ', options.RUNNUMBER

# Now fire up the random number generator with that seed
#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
    seed = 1234567,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.context['I3RandomService'] = randomService

def BasicHitFilter(frame):
    hits = 0
    if frame.Has(photon_series):
       hits = len(frame.Get(photon_series))
    if hits>0:
       return True
    else:
       return False


### START ###

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )


#tray.AddModule(ModMCTree, "modmctree", mctree="I3MCTree",
#	addhadrons=True
#	)

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

icemodel_path =  options.ICEMODEL
print 'Ice model ', icemodel_path
print "DOM efficiency: ", options.EFFICIENCY
# Only the photons are made. Still have to convert them to hits!
print "Setting cross energy: " , float(options.CROSSENERGY), "GeV"
#tray.AddSegment(clsim_hybrid.I3CLSimMakePhotons, 'goCLSIM',
print "Using CPUs ", CPU
print "Using GPUs ", options.GPU

gcd_file = dataio.I3File(options.GCDFILE)

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseCPUs=CPU,
                UseGPUs=options.GPU,
#		UseOnlyDeviceNumber=[0],
#                OpenCLDeviceList=[0],
                MCTreeName="I3MCTree",
                OutputMCTreeName="I3MCTree_clsim",
                FlasherInfoVectName=None,
                MMCTrackListName=None,
                PhotonSeriesName=photon_series,
                ParallelEvents=1000,
                RandomService=randomService,
                IceModelLocation=icemodel_path,
                #UnWeightedPhotons=True, #turn off optimizations
                UseGeant4=False,
                CrossoverEnergyEM=0.1,
		CrossoverEnergyHadron=float(options.CROSSENERGY),
                StopDetectedPhotons=True,
#                UseHoleIceParameterization=False, # Apply it when making hits!
               # HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.flasher_p1_0.30_p2_-1"),
                DoNotParallelize=False,
                DOMOversizeFactor=5.,
                UnshadowedFraction=options.EFFICIENCY, #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile=gcd_file,
                ExtraArgumentsToI3CLSimModule={
                    #"UseHardcodedDeepCoreSubdetector":True, #may save some GPU memory
                    #"EnableDoubleBuffering":True,
                    "DoublePrecision":False, #will impact performance if true
                    "StatisticsName":"clsim_stats",
                    "IgnoreDOMIDs":[],
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
