#!/usr/bin/env python

"""
Test Icetray code to run the noise generator
"""

import numpy as np
import argparse

import Salt

from icecube import icetray, dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from I3Tray import *
from icecube.dataclasses import ModuleKey


parser = argparse.ArgumentParser(description="Generates files of K40 noise hits")
parser.add_argument("-i", "--inFile", dest="inFile", help="path to input file location")
parser.add_argument(
    "-o", "--outFile", dest="outFile", help="path to output file location"
)
args = parser.parse_args()

tray = I3Tray()
icetray.logging.set_level(icetray.I3LogLevel.LOG_WARN)

tray.AddModule("I3Reader", "reader", FilenameList=[args.inFile])

tray.AddModule(
    Salt.K40Noise,
    "generatingNoise",
    charFile="/data/p-one/jstacho/noise/correlatedNoise/noiseGenerator/noiseModule/poneoffline/NoiseCharacterization.pkl",  # Path to noise characterization file
    inputTreeName="I3Photon",  # Name of the tree in the frames that determines how many frames are generated
    startPadding=1000,
    endPadding=10000,
    skipSingles=True,
    pmtHitsOnly=True,
    skipOneFold=True,
)

tray.AddModule(
    "I3Writer",
    "writer",
    Filename=args.outFile,
    Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
)

tray.AddModule("TrashCan", "adios")
tray.Execute()
tray.Finish()
