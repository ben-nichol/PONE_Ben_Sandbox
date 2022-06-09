#!/usr/bin/env python

# system imports
from optparse import OptionParser

# icetray imports
from I3Tray import I3Tray
from icetray import OMKey, I3Units, I3Frame
from icecube import dataio, phys_services, clsim

# pone imports
from Utilities.DOMUtility import DOMProperties
from Flasher.Isotropic import GenerateIsotropic
import WaterOpticalModel.MakePoneMediumPropertiesConservativeExtendedRange as Medium


# use parser options to setup simulation in prompt
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test-data/test-pocam.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12344,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default='',
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-fk", "--flasherkey", type="str", default='1-1',
                  dest="FLASHERKEY", help="Flasher position GCD string-om index")
parser.add_option("-ph", "--numphotons", type="int", default=1e5,
                  dest="NUMPHOTONS", help="The number of photons per flash")
parser.add_option("-w", "--fwhm", type="float", default=5,
                  dest="PULSEFWHM", help="Pulse FWHM  in nanoseconds")
parser.add_option("-oz", "--oversize", type="float", default=1,
                  dest="OVERSIZE", help="OM oversizing factor")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

# geometry
geometry = dataio.I3File(options.GCDFILE)
gframe = geometry.pop_frame()  
geo = gframe["I3Geometry"]

# flasher
fkey = [int(i) for i in options.FLASHERKEY.split('-')]
flasher_key = OMKey(fkey[0], fkey[1], 1)
flasher_position = geo.omgeo[flasher_key].position
flasher_photons = options.NUMPHOTONS
flasher_width = options.PULSEFWHM * I3Units.ns
flasher_pulse_type = clsim.I3CLSimFlasherPulse.FlasherPulseType.LED405nm

# dom properties instance and derived values
dom_properties = DOMProperties()
wl_acceptance_max = dom_properties.GetMaxAngularAcceptance() * 1.05
wl_acceptance = dom_properties.GetCLSimQETable(factor=wl_acceptance_max)
dom_radius = (17.0*2.54*0.01/2.0) * I3Units.m
dom_oversize = options.OVERSIZE

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
               Stream=I3Frame.DAQ)

# add event header
tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2022,
               DAQTime=7968509615844458,
               RunNumber=options.RUNNUMBER,
               EventID=1,
               IncrementEventID=True)

# add fake isotropic flasher similar to POCAM
tray.AddModule(GenerateIsotropic.GenerateIsotropic,
               SeriesFrameKey="FlasherPulseSeries",
               PhotonPosition=flasher_position,
               NumberOfPhotons=flasher_photons,
               PulseWidth=flasher_width,
               Seed=options.SEED,
      	       Isotropy=True,
               FlasherPulseType=flasher_pulse_type)

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
    FlasherPulseSeriesName="FlasherPulseSeries",
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
