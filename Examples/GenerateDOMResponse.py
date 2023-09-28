from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from I3Tray import *
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
import argparse
import DOM.PONEDOMLauncher 
#import PONEDOMLauncher 

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-n', '--runNum',  dest = 'runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-g', '--gcdFile', dest = 'GCDFile', help = "the GCD File used in the simulation")
parser.add_argument('-i', '--infile', dest = 'infile', help= 'input file and path')
parser.add_argument('-o', '--outfile', dest = 'outfile', help= 'output file and path')
args = parser.parse_args()

tray = I3Tray()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.GCDFile, args.infile]
            )

tray.AddModule(MCPEConverter.makeHits, "generatingHits",
		       AngularAcceptance = 0.811397114255,
               GCDFile = args.GCDFile,
               PropagatedPhotons = "I3Photons",
               OutputMCTreeName = "MCPESeriesMap",
               DOMOversizeFactor = 1.0,
               DefaultRelativeDOMEfficiency = 1.0)

tray.AddModule("I3Writer","writer",
               #SkipKeys=SkipKeys,
               Filename = args.outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")
tray.Execute()
tray.Finish()
