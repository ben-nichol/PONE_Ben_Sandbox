#!/usr/bin/env python

# import required icecube-related stuff
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Units
from I3Tray import I3Tray
# command line options required to configure the simulation
import argparse
from os.path import expandvars
import os, sys
from icecube import phys_services
#from icecube.simprod.modules import Corsika5ComponentGenerator
from segments.GenerateCosmicRayMuons import GenerateSingleMuons, GenerateNaturalRateMuons
from segments import GenerateCosmicRayMuons, PropagateMuons
#from icecube.simprod import segments
from Utilities.GeoUtility import get_geo_from_gcd

def printfunc(frame, message = 'test'):
	print(message)
	return True

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str,default="./test_output_muonprop.i3",help="")
parser.add_argument("-r", "--run",type=int,default=0,help="")                                                       
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Readin GCD file")
parser.add_argument("-n", "--nevents",type=int,default=1000,help="Number of events to run.")
parser.add_argument("-i", "--infile",type = str,default="./test_output.i3",help="")

args = parser.parse_args()

tray = I3Tray()

radius, height = get_geo_from_gcd(args.gcdfile)

randomService = phys_services.I3SPRNGRandomService(
                                  seed = args.run*args.run,
                              nstreams = 100000000,
                              streamnum = args.run)

tray.context['I3RandomService'] = randomService

_kwargs = {"PROPOSAL_config_file":os.getenv('PONESRCDIR')+"/configs/PROPOSAL_config.json"}

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile,args.infile]
            )

tray.Add(PropagateMuons, 'ParticlePropagators',
                                  RandomService=randomService,
                                  SaveState=True,
                                  InputMCTreeName="I3MCTree_preMuonProp",
                                  OutputMCTreeName="I3MCTree",
                                  PROPOSAL_config_file=os.getenv('PONESRCDIR')+"/configs/PROPOSAL_config.json")
#tray.Add(printfunc,"print4",message = 'print4')
tray.AddModule('I3Writer',   
                'writer',
                Streams=[icetray.I3Frame.Stream('S'),
                icetray.I3Frame.TrayInfo,
                icetray.I3Frame.DAQ],
                filename=args.outfile)
#tray.Add(printfunc,"print5",message = 'print5')
#tray.AddModule('TrashCan','YesWeCan')

tray.Execute()
tray.Finish()
