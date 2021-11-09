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
from segments.GenerateCosmicRayMuons import GenerateSingleMuons
from segments import GenerateCosmicRayMuons
#from icecube.simprod import segments

def printfunc(frame, message = 'test'):
	print(message)
	return True

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str,default="./test_output.root",help="")
parser.add_argument("-r", "--run",type=int,default=0,help="")                                                       
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Readin GCD file")
parser.add_argument("-n", "--nevents",type=int,default=10,help="Number of events to run.")

args = parser.parse_args()

tray = I3Tray()
tray.AddModule('I3InfiniteSource',Prefix=args.gcdfile)

randomService = phys_services.I3SPRNGRandomService(
                                  seed = args.run*args.run,
                              nstreams = 100000000,
                              streamnum = args.run)

tray.context['I3RandomService'] = randomService
#tray.Add(printfunc,"print0",message = 'print0')
#tray.AddModule('I3Reader', 'reader',  filenamelist = [args.gcdfile])
#tray.Add(printfunc,"print1",message = 'print1')
tray.Add("I3MCEventHeaderGenerator",
               EventID=1,
               IncrementEventID=True)
#tray.Add(printfunc,"print2",message = 'print2')
#tray.AddSegment(GenerateSingleMuons,"makeMuons",
#                GCDFile=args.gcdfile
#               )

#tray.AddSegment(Corsika5ComponentGenerator,"makeCRs",
#                nshowers=100000,
#                procnum=args.run,
#                nproc=1000,
#                seed=args.run*args.run,
#                gcdfile=args.gcdfile,
#                outputfile="pone_corsika_"+str(args.run)+".i3.zst",
#                RunCorsika=True,
#                sumaryfile="summary1.json",
#                pnorm=[5,2.25,1.1,1.2,1],
#                pgam=[2.65,2.6,2.6,2.6,2.6],
#                CORSIKAseed=args.run,
#                eprimarymax=100000,
#                eprimarymin=600,
#                OverSampling=1,
#                corsikaVersion="76900g",
#                CutoffType="EnergyPerParticle",
                #RepoURL="http://prod-exe.icecube.wisc.edu/",
#                UsePipe=True,
#                compress=True,
#                HistogramFilename="hist.pkl",
#                EnableHistogram=True
#               )    

tray.AddSegment(GenerateCosmicRayMuons,"CosmicRayMuons",
		mctree_name='I3MCTree_preMuonProp',
		num_events=args.nevents,
		flux_model='Hoerandel5_atmod12_SIBYLL',
		gamma_index=2.0,
		energy_offset=700.0,
		energy_min=1000.0,
		energy_max=1000000.0,
		cylinder_length=1080.0,
		cylinder_radius=240.0,
		cylinder_x=0.0,
		cylinder_y=0.0,
		cylinder_z=0.0,
#		inner_cylinder_length=500.0,
#		inner_cylinder_radius=150.0,
#		inner_cylindedr_x=46.3,
#		inner_cylinder_y=-34.9,
#		inner_cylinder_z=-300.0,
#		use_inner_cylinder=False
                )

tray.Add(segments.PropagateMuons, 'ParticlePropagators',
                                  RandomService=randomService,
                                  SaveState=True,
                                  InputMCTreeName="I3MCTree_preMuonProp",
                                  OutputMCTreeName="I3MCTree")
#tray.Add(printfunc,"print4",message = 'print4')
tray.AddModule('I3Writer',   
                'writer',
                Streams=[icetray.I3Frame.Stream('S'),
                icetray.I3Frame.TrayInfo,
                icetray.I3Frame.DAQ],
                filename=args.outfile+"Corsika_lowE"+str(args.run)+".i3.gz")
#tray.Add(printfunc,"print5",message = 'print5')
#tray.AddModule('TrashCan','YesWeCan')

tray.Execute()
tray.Finish()
