#!/bin/sh bash
#METAPROJECT combo/V00-00-04

from os.path import expandvars
import os, sys, random
from I3Tray import *                                                            
import random                                                                   
from icecube import icetray, dataclasses, dataio, simclasses                    
from icecube import phys_services, sim_services           
import argparse  
from Reconstruction.Linefit.LineFitReco import LineFitReco
from Reconstruction.Track.TrackReco import TrackReco
from PulseCleaning.SignificantHitPulseCleaning import SignificantHitPulseCleaning
from PulseCleaning.CausalHits import CausalPulseCleaning
from Reconstruction.Cascade.CascadeReco import CascadeReco
from Reconstruction.TrackLimits.StartStop import StartStopFit

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="./test_output.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type=int, default="1", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",type=int,default=1, help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Read in GCD file")
parser.add_argument("-t", "--pulsesep",default=0.000001,help="Time needed to separate two pulses. Assume that this is 3.5*sample time.")
parser.add_argument("-e", "--ext",default=".zst",help="compression extension")

args = parser.parse_args()
photon_series = "I3Photons"
tray = I3Tray()

#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                                                   seed = 1234567,
                                                   nstreams = 10000,
                                                   streamnum = args.runnumber
                                                   )

tray.context['I3RandomService'] = randomService

infile = args.infile +str(args.runnumber)+".i3"+args.ext
outfile = args.outfile +str(args.runnumber)+".i3"+args.ext

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile, infile]
            )

print(args.gcdfile)
gcd_file = dataio.I3File(args.gcdfile)

#This pulse cleaning is promissing but still experimental. 
tray.AddModule(CausalPulseCleaning,"CausalHit",
              GCDFile=gcd_file,
              inputseries = "PMTResponse",
              output = "CausalHits"
              )

#Linefit for tracks
tray.AddModule(LineFitReco, "LineFit",
              inputseries = "CausalHits",
              output = "linefit"
              )

#Track reconstruction
tray.AddModule(TrackReco,"likelihoodreco",
               pulseseries = "CausalHits",
               seedtrack = "linefit",
               output = "llhfit",
              )

tray.AddModule(StartStopFit,"StartStop",
                pulseseries = "CausalHits",
                seedtrack = "llhfit",
                output = "startstop"
              )

tray.AddModule(CascadeReco,"NuTauReconstruction",
              pulseseries = "CausalHits",
              output = "Cascade",
              )

tray.AddModule("I3Writer","writer",
               Filename = outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
