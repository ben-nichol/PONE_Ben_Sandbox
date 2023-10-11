#**************************************
# P-ONE example to propagate muons from the generator
# through the P-ONE detector medium
#**************************************

# import required icecube-related stuff
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Units, I3Tray
# command line options required to configure the simulation
import argparse
import os, sys
from os.path import expandvars
from icecube import phys_services
from Utilities.GeoUtility import get_geo_from_gcd

#Makes use of a predefined "segement" to do the heavy lifting
from segments import PropagateMuons

#Read in user options, or set to default
parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str,default="dataio/muonprop.i3",help="")
parser.add_argument("-r", "--run",type=int,default=0,help="")                                                       
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_5String.i3.gz", help="Readin GCD file")
parser.add_argument("-n", "--nevents",type=int,default=1000,help="Number of events to run.")
parser.add_argument("-i", "--infile",type = str,default="dataio/data_output.i3",help="")

args = parser.parse_args()

#Initial setup
tray = I3Tray()

radius, height = get_geo_from_gcd(args.gcdfile)

randomService = phys_services.I3SPRNGRandomService(
                                  seed = args.run*args.run,
                              nstreams = 100000000,
                              streamnum = args.run)

tray.context['I3RandomService'] = randomService


tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile,args.infile]
            )

tray.Add(PropagateMuons, 'ParticlePropagators',
                                  RandomService=randomService,
                                  SaveState=True,
                                  InputMCTreeName="I3MCTree",
                                  OutputMCTreeName="I3MCTree",
                                  PROPOSAL_config_file=os.getenv('PONESRCDIR')+"/configs/PROPOSAL_config.json")

tray.AddModule('I3Writer',   
                'writer',
                Streams=[icetray.I3Frame.Stream('S'),
                icetray.I3Frame.TrayInfo,
                icetray.I3Frame.DAQ],
                filename=args.outfile)

tray.Execute()
tray.Finish()
