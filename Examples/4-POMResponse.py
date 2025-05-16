#!/bin/sh

import argparse
import os, sys, random
from os.path import expandvars

from DOM.DOMAcceptance import DOMAcceptance
from DOM.PONEDOMLauncher import DOMSimulation

from Trigger.DOMTrigger import DOMTrigger

from NoiseGenerators.DarkNoise import DarkNoise
from NoiseGenerators.K40Noise import K40Noise

from icecube.icetray import I3Tray
from icecube import phys_services, sim_services
from icecube import icetray, dataclasses, dataio, simclasses

from PulseCleaning.CausalHits import CausalPulseCleaning
from Trigger.DOMTrigger import DOMTrigger
from Trigger.DetectorTrigger import DetectorTrigger

parser = argparse.ArgumentParser()
parser.add_argument(
    "-o",
    "--outfile",
    type=str,
    default="dataio/POM_response.i3",
    help="Write output to OUTFILE (.i3{.gz} format)",
)
parser.add_argument(
    "-i",
    "--infile",
    type=str,
    default="dataio/clsim.i3",
    help="Read input from INFILE (.i3{.gz} format)",
)
parser.add_argument(
    "-r",
    "--runnumber",
    type=int,
    default="1",
    help="The run/dataset number for this simulation, is used as seed for random generator",
)
parser.add_argument(
    "-l",
    "--filenr",
    type=int,
    default=1,
    help="File number, stream of I3SPRNGRandomService",
)
parser.add_argument(
    "-g",
    "--gcdfile",
    default=os.getenv("PONESRCDIR") + "/GCD/PONE_5String.i3.gz",
    help="Read in GCD file",
)
parser.add_argument(
    "-t",
    "--pulsesep",
    default=0.2,
    help="Time needed to separate two pulses. Assume that this is 3.5*sample time.",
)
parser.add_argument("-e", "--ext", default=".gz", help="compression extension")
parser.add_argument(
    "-s",
    "--dropstrings",
    nargs="+",
    default=[],
    help="Strings to exclude from geometry",
)
parser.add_argument(
    "-n", "--nDOMs", type=int, default=1, help="Number of DOMs for detector trigger"
)
parser.add_argument(
    "-f",
    "--LICconfig",
    type=str,
    default="",
    help="Path to the LIC configuration file for Lepton Injection events.",
)
parser.add_argument(
    "-c",
    "--crossdir",
    default=os.getenv("PONESRCDIR") + "/CrossSectionModels/csms_differential_v1.0",
    help="path to cross section models",
)

tray = I3Tray()

args = parser.parse_args()
photon_series = "I3Photons"
tray = I3Tray()

dropstrings = []
for string in args.dropstrings:
    dropstrings.append(int(string))

# from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
    seed=1234567, nstreams=10000, streamnum=args.runnumber
)

tray.context["I3RandomService"] = randomService

infile = args.infile
outfile = args.outfile

tray.AddModule("I3Reader", "reader", FilenameList=[args.gcdfile, infile])

tray.AddModule(DOMAcceptance,
               'DOMAcceptance',
               input_map      = photon_series,
               output_map     = 'Accepted_PulseMap',
               random_service = randomService
               )


tray.AddModule(DarkNoise,
               'AddDarkNoise',
               input_map      = 'Accepted_PulseMap',
               output_map     = 'Noise_Dark',
               random_service = randomService,
               drop_oms       = [2, 3],
               gcd_file       = os.getenv('PONESRCDIR') + '/GCD/one-om-gcd-origin.i3.gz'
               )


tray.AddModule(K40Noise,
               'AddK40Noise',
               input_map             = 'Accepted_PulseMap',
               output_map            = 'Noise_K40',
               characterization_file = os.getenv('PONESRCDIR') + '/NoiseGenerators/k40-characterization.pkl',
               random_service        = randomService,
               drop_oms              = [2, 3],
               gcd_file              = os.getenv('PONESRCDIR') + '/GCD/one-om-gcd-origin.i3.gz'
               )


tray.AddModule(DOMSimulation,
               'DOMLauncher',
               input_map      = 'Accepted_PulseMap',
               output_map     = 'PMT_Response',
               random_service = randomService,
               min_time_sep   = args.pulsesep,
               split_doms     = True,
               use_dark       = True,
               dark_map       = 'Noise_Dark',
               use_k40        = True,
               k40_map        = 'Noise_K40'
              )

tray.AddModule(
    DOMTrigger,
    "DOMTrigger",
    inputmap="PMTResponse",
)

tray.AddModule(
    DetectorTrigger,
    "PONE_Trigger",
    output="_3PMT_2DOM",
    DOMPMTCoinc=3,
    FullDetectorCoincidenceN=args.nDOMs,
    CutOnTrigger=True,
    EventLength=10000,
    TriggerTime=2000,
    PulseSeriesIn="PMTResponse",
    PulseSeriesOut="EventPulseSeries",
)

tray.AddModule(
    "I3Writer",
    "writer",
    # SkipKeys = ["I3Photons","I3Photons_PMTResponse","TimeShiftedMCPEMap"],
    Filename=outfile,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
)

tray.AddModule("TrashCan", "adios")

tray.Execute()
tray.Finish()
