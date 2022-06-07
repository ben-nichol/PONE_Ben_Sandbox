#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

# use parser options to setup simulation in prompt
usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test-flashers.i3",
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

# system imports
import os
import sys
import math
import numpy

# icetray imports
from I3Tray import *
from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services

# pone imports
import WaterOpticalModel.MakePoneMediumPropertiesConservativeExtendedRange as Medium
from Utilities.DOMUtility import DOMProperties

# dom properties instance and derived values
dom_properties = DOMProperties()
wl_acceptance_max = dom_properties.GetMaxAngularAcceptance() * 1.05
wl_acceptance = dom_properties.GetCLSimQETable(factor=wl_acceptance_max)
dom_radius = (17.0*2.54*0.01/2.0)*icetray.I3Units.m
dom_oversize = 1.0 # 17./13.

# optical medium
optical_medium = Medium.MakePoneMediumProperties()

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

# start icecube tray
tray = I3Tray()

# add geometry and daq stream
tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

# add event header
tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2012,
               DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

# add fake IceCube LED flashers
# flasher mask = hhhhhhpppppp (h=horizontal, p=polar LEDs)
# 1=enable, 0=disable
# brightness [0, 127], low to high
# width [0, 127], low to high
tray.AddModule(clsim.FakeFlasherInfoGenerator, "FakeFlasherInfoGenerator",
               FlashingDOM = icetray.OMKey(5, 10, 1), # flasher position
               FlasherTime = 0.*I3Units.ns, # start time of flash
               FlasherMask = 0b100000000000, # only the 1st horizontal LED 
               FlasherBrightness = 60, # brightness
               FlasherWidth = 10)      # width, rectangular pulse

# start photon propagation with CLsim
tray.AddSegment(clsim.I3CLSimMakePhotons, "goCLSIM",
    UseGPUs=True,
    UseCPUs=False,
    UseGeant4=False,
    UseI3PropagatorService=False,
    RandomService=randomService,
    DoNotParallelize=False,
    UnweightedPhotons=True,
    StopDetectedPhotons=True,
    PhotonSeriesName = 'I3Photons',
    MCPESeriesName='',
    FlasherInfoVectName="I3FlasherInfo",
    GCDFile=options.GCDFILE,
    IceModelLocation=optical_medium,
    WavelengthAcceptance=wl_acceptance,
    DOMRadius=dom_radius,
    DOMOversizeFactor=dom_oversize,
)

# write propagated photons to file
tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)

# execute
tray.Execute(options.NUMEVENTS + 3)
