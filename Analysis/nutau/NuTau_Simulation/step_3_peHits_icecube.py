import sys
sys.path.insert(0,'/home/users/dhilu/P_ONE_dvirhilu/src')
sys.path.insert(0, '/home/users/jguthrie/oscnext/oscnext_scripts')

from optparse import OptionParser
from os.path import expandvars
import os, random

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
#parser.add_option("-r", "--runnumber", type="string", default="1",
#                  dest="RUNNUMBER", help="The run number for this simulation, is used as seed for random generator")
#parser.add_option("-f", "--filenr",type="string",default="1",
#                  dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_option("-g", "--gcdfile", default=os.getenv('GCDfile'),
		          dest="GCDFILE", help="Read in GCD file")
parser.add_option("-e","--efficiency", type="float",default=1.0,
                  dest="EFFICIENCY",help="DOM Efficiency ... the same as UnshadowedFraction")
#parser.add_option("-n","--noise", default="vuvuzela",
#                  dest="NOISE",help="Noise model (vuvuzela/poisson)")
#parser.add_option("-l", "--holeice",  default = "as.flasher_p1_0.30_p2_-1",
#                  dest="HOLEICE",
#                  help="Pick the hole ice parameterization, corresponds to a file name in $I3_SRC/ice-models/resources/models/angsens/")
parser.add_option("-m","--icemodel", default="spice_3.2.1",
                  dest="ICEMODEL",help="Ice model and path")

(options,args) = parser.parse_args()

icemodel_path = options.ICEMODEL

from I3Tray import *
import random

from icecube import icetray, dataclasses, dataio, simclasses #, recclasses
from icecube import phys_services, sim_services, DOMLauncher, DomTools, clsim, trigger_sim

from RemoveLatePhotons_V5 import RemoveLatePhotons

def BasicHitFilter(frame):
    hits = 0
    if frame.Has("MCPESeriesMap"):
       hits = len(frame.Get("MCPESeriesMap"))
    if hits>0:
#        print "has photons"
        return True
    else:
#       print "does NOT photons"
        return False

def BasicDOMFilter(frame):
    if frame.Has("InIceRawData"):
        if len(frame['InIceRawData']) > 0:
            return True
        else:
            return False
    else:
       return False

tray = I3Tray()

randomService = phys_services.I3GSLRandomService(123456)
tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom")(
    ("Seed", 123456)
    )

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )

#tray.AddModule(RemoveLatePhotons, "RemovePhotons",
               #InputPhotonSeries = "I3Photons",
               #TimeLimit = 1E5) #nanoseconds

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")
gcd_file = dataio.I3File(options.GCDFILE)

tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons, "makeHitsFromPhotons",
#                MCTreeName="I3MCTree_clsim",
#                PhotonSeriesName="UnweightedPhotons2",
                PhotonSeriesName="I3Photons",
                MCPESeriesName="MCPESeriesMap",
                RandomService=randomService,
                DOMOversizeFactor=5.0,
                UnshadowedFraction=options.EFFICIENCY,
                IceModelLocation = icemodel_path,
#               UseHoleIceParameterization=holeice
#               HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/%s"%options.HOLEICE),
                GCDFile=gcd_file
                ) #generates PEs

mcpe_to_pmt = "MCPESeriesMap"
mcpeout = mcpe_to_pmt

#tray.AddModule("PMTResponseSimulator","rosencrantz",
    #Input=mcpeout,
    #Output=mcpeout + "_weighted",
    #MergeHits=True,
    #RandomServiceName=RandomService,
    #) -----> Commented Out, generates pulses

#tray.AddModule("DOMLauncher", "guildenstern",
    #Input= mcpeout + "_weighted",
    #Output="InIceRawData_unclean",
    #UseTabulatedPT=True,
    #)

#tray.AddModule("I3DOMLaunchCleaning","launchcleaning")(
       #("InIceInput","InIceRawData_unclean"),
       #("InIceOutput","InIceRawData"),
       #("FirstLaunchCleaning",False),
#       ("CleanedKeys",BadDoms)
       #) -----> Commented Out:

tray.AddModule(BasicHitFilter, 'FilterNullInIce', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

tray.AddModule('Delete', 'delete_triggerHierarchy',
               Keys = ['I3TriggerHierarchy', 'TimeShift', 'CleanIceTopRawData'])

#tray.AddSegment(trigger_sim.TriggerSim, 'trig',
                #gcd_file = gcd_file,
 #               time_shift_args = time_shift_args,
#added in run_id
                #run_id=1)

# Not skipping these keys for now (check what gets dropped in the L2)
skipkeys = ["MCPMTResponseMap",
            "MCTimeIncEventID",
            "I3MCTree_clsim"
            "I3Photons"
            "clsim_stats",
            "InIceRawData_unclean",
            ]

tray.AddModule("I3Writer","writer",
               #SkipKeys=skipkeys, All of these get thrown out by the L2 anyways ... keep them?
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
