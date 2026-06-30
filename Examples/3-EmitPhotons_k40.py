#!/usr/bin/env python3


from icecube import (
    icetray,
    dataclasses,
    dataio,
    phys_services,
    sim_services,
    simclasses,
    clsim,
)
from icecube.icetray import I3Tray, I3Units
from argparse import ArgumentParser
from os.path import expandvars
from Utilities.DOMUtility import DOMProperties
import os
import WaterOpticalModel.MakePoneMediumPropertiesConservative as Medium


usage = "usage: %prog [options] inputfile"
parser = ArgumentParser(usage)
parser.add_argument(
    "-i",
    "--infile",
#    default="dataio/prop.i3",
    default="k40test.i3.gz",
    help="Write output to OUTFILE (.i3{.gz} format)",
)
parser.add_argument(
    "-o",
    "--outfile",
    default="dataio/clsim_k40.i3",
    help="Write output to OUTFILE (.i3{.gz} format)",
)
parser.add_argument(
    "-g", "--gcdfile", default=os.getenv("PONESRCDIR") + "/GCD/PONE_5String.i3.gz"
)
parser.add_argument(
    "--icemodel",
    default=expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_lea"),
    dest="ICEMODEL",
    help="A clsim ice model file/directory (ice models *will* affect performance metrics, always compare using the same model!)",
)
parser.add_argument(
    "--use-cpu",
    action="store_true",
    default=False,
    dest="USECPU",
    help="simulate using CPU instead of GPU",
)

# parse cmd line args, bail out if anything is not understood
options = parser.parse_args()

dom_properties = DOMProperties()

tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed=int(11234), nstreams=int(4e7), streamnum=int(452)
)

tray.context["I3RandomService"] = randomService

outfile = options.outfile
infile = options.infile

tray.AddModule("I3Reader", "reader", FilenameList=[infile])

MCTreeName = "I3DecayProducts"
photonSeriesName = "I3Photons"

# NEW THINGS TO MAKE WORK FOR 1.14
WavelengthAcceptance = dom_properties.GetCLSimQETable(factor=dom_properties.GetMaxAngularAcceptance() * 1.05)
GCDFile  = options.gcdfile
infile=dataio.I3File(GCDFile)
frames = []
for frame in infile:
	frames.append(frame)
calib_frame = frames[1]
rde = dict()
for k, domcal in dataclasses.I3Calibration.from_frame(calib_frame).dom_cal.items():
	rde[k] = domcal.relative_dom_eff
domAcceptance = simclasses.I3CLSimFunctionMap()
for k in rde.keys():
	domAcceptance[k] = WavelengthAcceptance
# END OF NEW THINGS


kwargs = {}

tray.AddSegment(
    clsim.I3CLSimMakeHits,
    "makeCLSimHits",
    GCDFile=options.gcdfile,
    DOMRadius=0.21590 * icetray.I3Units.m,  # POM dimension, default is icecube specific
    IceModelLocation=Medium.MakePoneMediumProperties(),
    MCPESeriesName="",
    PhotonSeriesName=photonSeriesName,
    MCTreeName=MCTreeName,
    RandomService=randomService,
    DOMEfficiency=0.95,
    UseGPUs=not options.USECPU,
    UseCPUs=options.USECPU,
    StopDetectedPhotons=True,
    HoleIceParameterization=os.getenv("PONESRCDIR") + "/data/as.full",
    DoNotParallelize=False,
    UnshadowedFraction=1.0,  # normal in IC79 and older CLSim versions was 0.9, now it is 1.0
    UseI3PropagatorService=False,
    UnWeightedPhotons=True,
    WavelengthAcceptance=domAcceptance,
)

tray.AddModule(
    "I3Writer",
    "writer",
    #               SkipKeys=SkipKeys,
    Filename=outfile,
    Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation, icetray.I3Frame.DAQ],
)


tray.AddModule("TrashCan", "adios")

tray.Execute()
tray.Finish()
