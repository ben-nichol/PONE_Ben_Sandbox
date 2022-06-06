#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test-flashes.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12344,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default=expandvars("$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services

import WaterOpticalModel.MakePoneMediumPropertiesConservativeExtendedRange as Medium
from Utilities.DOMUtility import DOMProperties

import math
import numpy

dom_properties = DOMProperties()

tray = I3Tray()

# a random number generator
try:
    randomService = phys_services.I3SPRNGRandomService(
        seed = options.SEED,
        nstreams = 10000,
        streamnum = options.RUNNUMBER)
except AttributeError:
    randomService = phys_services.I3GSLRandomService(
        seed = options.SEED*10000 + options.RUNNUMBER,
    )

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2012,
               DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

tray.AddModule(clsim.FakeFlasherInfoGenerator, "FakeFlasherInfoGenerator",
               FlashingDOM = icetray.OMKey(5, 10, 1),
               FlasherTime = 0.*I3Units.ns,
               FlasherMask = 0b111111000000, # only the 6 horizontal LEDs 
               FlasherBrightness = 127, # full brightness
               FlasherWidth = 127)      # full width

#tray.AddSegment(clsim.I3CLSimMakePhotons, "goCLSIM",
#    UseGPUs=True,
#    UseCPUs=False,
#    UseGeant4=False,
#    UseI3PropagatorService=False,
#    RandomService = randomService,
#    DoNotParallelize=False,
#    UnweightedPhotons=False,
#    StopDetectedPhotons=True,
#    PhotonSeriesName = 'I3Photons',
#    MCPESeriesName='',
#    FlasherInfoVectName="I3FlasherInfo",
#    GCDFile = options.GCDFILE,
#    IceModelLocation=Medium.MakePoneMediumProperties(),
#    WavelengthAcceptance = dom_properties.GetCLSimQETable( factor=dom_properties.GetMaxAngularAcceptance()*1.05 ),
#    DOMRadius = (17.0*2.54*0.01/2.0)*icetray.I3Units.m,
#    DOMOversizeFactor=1.0, #(17./13.),
#)

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute(options.NUMEVENTS+3)
