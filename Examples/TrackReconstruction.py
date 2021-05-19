#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT combo/V00-00-04

from os.path import expandvars
import os, sys, random
from DOM.PONEDOMLauncher import SimpleDOMSimulation
from DOM.WaveformBuilder import WaveformBuilder
from PulseCleaning.TimeShifted import timeShift
from I3Tray import *                                                            
import random                                                                   
from icecube import icetray, dataclasses, dataio, simclasses                    
from icecube import phys_services, sim_services           
import argparse  
from Reconstruction.linefit.SimAnalysis import LineFitReco
from Reconstruction.llh.likelihoodreco_tables import likelihoodreco
from Reconstruction.nutau.NuTauReco_tables import nutaureco
from PulseCleaning.SignificantHitPulseCleaning import SignificantHitPulseCleaning
from Reconstruction.nutau.curveFit import curveFit
from Reconstruction.nutau.curveFit_tables import curveFit_tables 

# This script will perform a hybridCLSim propagation.
#
# NOTE: There is no bad_dom_cleaning!!!
#       This you still have to do after the propagation!!!

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="./test_output.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type=int, default="1", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",type=int,default=1, help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="Read in GCD file")

args = parser.parse_args()
photon_series = "I3Photons"
tray = I3Tray()

files_dir = args.infile
file_list_aux = os.listdir(files_dir)
file_list = [x for x in file_list_aux if ( '.i3.gz' in x)]

#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                                                   seed = 1234567,
                                                   nstreams = 10000,
                                                   streamnum = args.runnumber
                                                   )

tray.context['I3RandomService'] = randomService

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile, args.infile+file_list[args.runnumber]]
            )

gcd_file = dataio.I3File(args.gcdfile)

#Take the MC photons and shifts them so that all events have a similar time structure.
tray.AddModule(timeShift,"MCtimeShift",
              MergedMCPETreeName = photon_series,
              TimeShiftedMCPE = "TimeShiftedMCPEMap",
              MinTime = 7200
              )

#Takes the MC photons and simulates to characteristics of the DOMs.
tray.AddModule(SimpleDOMSimulation, 'DOMLauncher',
               GCDFile=gcd_file,
               inputmap = "TimeShiftedMCPEMap",
               outputmap = "I3Photons_PMTResponse",
               RandomService = randomService
              )

#This is not typically needed but is for building raw waveforms for testing DOMLauncher. 
tray.AddModule(WaveformBuilder,'waveformbuilder',
              inputmap = "I3Photons_PMTResponse_MCpulses",
              outputmap = "waveform_pulses")

#This is a first order pulse cleaning method that will just remove pulses with time residualts > 1000 ns. 
tray.AddModule(SignificantHitPulseCleaning,"SignificantHit",
              GCDFile=gcd_file,
              inputseries = "I3Photons_PMTResponse",
              output = "SignificanHits",
              window = 1000
              )

#First order track reconstruction. 
tray.AddModule(LineFitReco, "LineFit",
              inputseries = "SignificanHits",
              output = "linefit"
              )


#Likelihood based track reconstruction. 
tray.AddModule(likelihoodreco,"likelihoodreco",
               pulseseries = "SignificanHits",
               seedtrack = "linefit",
               output = "llhfit"
              ) 

#An incomplete cascade reconstruction. 
tray.AddModule(nutaureco,"NuTauReconstructin",
              pulseseries = "SignificanHits",
              output = "NuTau"
              )

tray.AddModule("I3Writer","writer",
               Filename = args.outfile+"/NuTauFit_"+file_list[args.runnumber],
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
