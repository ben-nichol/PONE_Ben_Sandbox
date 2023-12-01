#!/bin/sh 

from os.path import expandvars
import os, sys, random
from DOM.PONEDOMLauncher import SimpleDOMSimulation
#from DOM.DAQTimeErrorSim import TriggerPulseErrorSim
import random                                                                   
from icecube import icetray, dataclasses, dataio, simclasses                    
from icecube import phys_services, sim_services           
from icecube.icetray import I3Tray           
import argparse  
from PulseCleaning.CausalHits import CausalPulseCleaning
from Trigger.DOMTrigger import DOMTrigger
from Trigger.DetectorTrigger import DetectorTrigger
#from Weighting.LeptonWeighter import LeptonWeighter

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="dataio/POM_response.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-i", "--infile",type=str, default="dataio/photons.i3", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type=int, default="1", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",type=int,default=1, help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_5String.i3.gz", help="Read in GCD file")
parser.add_argument("-t", "--pulsesep",default=0.2,help="Time needed to separate two pulses. Assume that this is 3.5*sample time.")
parser.add_argument("-e", "--ext",default=".gz",help="compression extension")
parser.add_argument("-s", "--dropstrings",nargs="+",default=[],help="Strings to exclude from geometry")
parser.add_argument("-n", "--nDOMs",type=int,default=1, help="Number of DOMs for detector trigger")
parser.add_argument("-f", "--LICconfig",type=str,default="",help="Path to the LIC configuration file for Lepton Injection events.")
parser.add_argument("-c", "--crossdir",  default=os.getenv('PONESRCDIR')+"/CrossSectionModels/csms_differential_v1.0",    help='path to cross section models')

tray = I3Tray()

args = parser.parse_args()
photon_series = "I3Photons"
tray = I3Tray()

dropstrings = []
for string in args.dropstrings :
    dropstrings.append(int(string))

#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                                                   seed = 1234567,
                                                   nstreams = 10000,
                                                   streamnum = args.runnumber
                                                   )

tray.context['I3RandomService'] = randomService

infile = args.infile
outfile = args.outfile

#if os.path.isfile(outfile) :
#    if os.path.getmtime(infile) < os.path.getmtime(outfile) :
#        quit()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.gcdfile, infile]
            )

tray.AddModule(SimpleDOMSimulation, 'DOMLauncher',
               inputmap = photon_series,
               outputmap = "PMTResponse",
               RandomService = randomService,
               minTsep = args.pulsesep,
               SplitDoms = True,
               dropstrings = dropstrings,
               add_noise = True,
               #DOMAcceptanceFile='/home/users/ctung/data/numu0002/jobfiles/PMTAcceptance_13PMTConfig.txt',
              )

tray.AddModule(DOMTrigger,"DOMTrigger",
               inputmap = 'PMTResponse',
              )

tray.AddModule(DetectorTrigger,"PONE_Trigger",
               output="_3PMT_2DOM",
               DOMPMTCoinc = 3,
               FullDetectorCoincidenceN = args.nDOMs,
               CutOnTrigger = True,
               EventLength = 10000,
               TriggerTime = 2000,
               PulseSeriesIn = "PMTResponse",
               PulseSeriesOut = "EventPulseSeries"
              )

tray.AddModule("I3Writer","writer",
               #SkipKeys = ["I3Photons","I3Photons_PMTResponse","TimeShiftedMCPEMap"],
               Filename = outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
