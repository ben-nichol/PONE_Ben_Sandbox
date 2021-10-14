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
from Reconstruction.Linefit.LineFitReco import LineFitReco
from Reconstruction.Track.TrackReco import TrackReco
from PulseCleaning.SignificantHitPulseCleaning import SignificantHitPulseCleaning
from Reconstruction.Cascade.CascadeReco import CascadeReco
from Trigger.DOMTrigger import DOMTrigger
from Trigger.DetectorTrigger import DetectorTrigger
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
parser.add_argument("-t", "--pulsesep",default=0.000001,help="Time needed to separate two pulses. Assume that this is 3.5*sample time.")

args = parser.parse_args()
photon_series = "I3Photons"
tray = I3Tray()

files_dir = args.infile
file_list_aux = os.listdir(files_dir)
file_list = [x for x in file_list_aux if ( '.i3.gz' in x and 'PhotonProp' in x and 'Reco' not in x)]

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

print(args.gcdfile)
gcd_file = dataio.I3File(args.gcdfile)

tray.AddModule(timeShift,"MCtimeShift",
              MergedMCPETreeName = photon_series,
              TimeShiftedMCPE = "TimeShiftedMCPEMap",
              MinTime = 7200
              )

tray.AddModule(SimpleDOMSimulation, 'DOMLauncher',
               GCDFile=gcd_file,
               inputmap = "TimeShiftedMCPEMap",
               outputmap = "I3Photons_PMTResponse",
               RandomService = randomService,
               minTsep = args.pulsesep,
               LPprob = 0.00000001,
               DNprob = 0.000000000001,
               APprob = 0.000000000001,
               SplitDoms = True,
               DOMAcceptanceFile = "/home/users/tmcelroy/pone_offline/data/PMTAcceptance_13PMTConfig.txt",
               PMTQEFile = "/home/users/tmcelroy/pone_offline/data/PMTQE.txt",
               AcceptBaseValue = -1.0
              )

tray.AddModule(DOMTrigger,"DOMTrigger",
                GCDFile=gcd_file,
                inputmap = "I3Photons_PMTResponse",
              )

tray.AddModule(DetectorTrigger,"PONE_Trigger",
               GCDFile=gcd_file,
               output="_2PMT_4DOM",
               DOMPMTCoinc =2,
               FullDetectorCoincidenceN = 4,
               CutOnTrigger = True
              )

tray.AddModule(SignificantHitPulseCleaning,"SignificantHit",
              GCDFile=gcd_file,
              inputseries = "I3Photons_PMTResponse",
              output = "SignificanHits",
              window = 1000
              )

tray.AddModule(LineFitReco, "LineFit",
              inputseries = "SignificanHits",
              output = "linefit"
              )

tray.AddModule(TrackReco,"likelihoodreco",
               pulseseries = "SignificanHits",
               seedtrack = "linefit",
               output = "llhfit",
              )

tray.AddModule(CascadeReco,"NuTauReconstruction",
              pulseseries = "SignificanHits",
              output = "Cascade",
              )

tray.AddModule("I3Writer","writer",
               Filename = args.outfile+"/TrigReco_"+file_list[args.runnumber],
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
