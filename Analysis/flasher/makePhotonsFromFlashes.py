 #python

from optparse import OptionParser
from os.path import expandvars
import os, sys, random

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_photons.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12344,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-m","--icemodel", default="spice_3.2.1",
                  dest="ICEMODEL",help="Ice model (spice_mie, spice_lea, etc)")
parser.add_option("-t", action="store_true",  dest="GPU", default=True,help="Run on GPUs or CPUs")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

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

def BasicHitFilter(frame):
    hits = 0
    if frame.Has(photon_series):
       hits = len(frame.Get(photon_series))
    if hits>0:
       return True
    else:
       return False

#print 'CUDA devices: ', options.DEVICE
tray = I3Tray()
print 'Using RUNNUMBER: ', options.RUNNUMBER

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

#icemodel_path = '/home/users/akatil/software/V06-01-02/build/ice-models/resources/models/' + str(options.ICEMODEL)
icemodel_path = str(options.ICEMODEL)
print 'Ice model ', icemodel_path
#print "DOM efficiency: ", options.EFFICIENCY
# Only the photons are made. Still have to convert them to hits!
#print "Setting cross energy: " , float(options.CROSSENERGY), "GeV"
#tray.AddSegment(clsim_hybrid.I3CLSimMakePhotons, 'goCLSIM',
print "Using CPUs ", CPU
print "Using GPUs ", options.GPU

gcd_file = dataio.I3File(options.GCDFILE)
#mediumProperties = clsim.MakeAntaresMediumProperties()

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseCPUs=CPU,
                UseGPUs=options.GPU,
		        #UseOnlyDeviceNumber=[1],
                #OpenCLDeviceList=[0],
                #MCTreeName="I3MCTree",
                OutputMCTreeName="I3MCTree_clsim",
                #FlasherInfoVectName="I3FlasherInfo",
                FlasherPulseSeriesName="CustomFlashes",
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
                UnshadowedFraction=1.2, #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
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

#tray.AddModule(BasicHitFilter, 'FilterNullPhotons', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

SkipKeys = ["I3MCTree_bak"]

tray.AddModule("I3Writer","writer",
               SkipKeys=SkipKeys,
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
