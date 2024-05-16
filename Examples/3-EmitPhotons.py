#!/usr/bin/env python3


from icecube import icetray, dataclasses, dataio, phys_services, sim_services, simclasses, clsim
from icecube.icetray import I3Tray, I3Units
from argparse import ArgumentParser
from os.path import expandvars
from Utilities.DOMUtility import DOMProperties
import os
import WaterOpticalModel.MakePoneMediumPropertiesConservative as Medium


usage = "usage: %prog [options] inputfile"
parser = ArgumentParser(usage)
parser.add_argument("-i", "--infile",default="dataio/prop.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-o", "--outfile",default="./clsim.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-x", "--xmlfile", default=None,dest="JSONFILE", help="Write statistics to JSONFILE")
parser.add_argument("--oversize", default=1, type=float,dest="OVERSIZE", help="DOM oversize factor")
parser.add_argument("--unweighted-photons", action="store_true",help="Propagate all Cherenkov photons. This is ~13x slower than downsampling first.")
parser.add_argument("-g", "--gcdfile",default=os.getenv("PONESRCDIR") + "/GCD/PONE_5String.i3.gz")
#parser.add_argument("-g", "--gcdfile",default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withStdNoise.i3.gz")
parser.add_argument("--use-cpu",  action="store_true", default=False,dest="USECPU", help="simulate using CPU instead of GPU")
parser.add_argument("--double-buffering", default=False, action="store_true",help="Interleave kernel execution and i/o")
parser.add_argument("--icemodel", default=expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_lea"),dest="ICEMODEL", help="A clsim ice model file/directory (ice models *will* affect performance metrics, always compare using the same model!)")

# parse cmd line args, bail out if anything is not understood
options = parser.parse_args()
# icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_INFO)
#icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_WARN)

dom_properties = DOMProperties()

tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
                seed = int(11234),
               nstreams = int(4e7),
                streamnum = int(452))

tray.context['I3RandomService'] = randomService

outfile = options.outfile
infile = options.infile

tray.AddModule('I3Reader', 'reader',
            FilenameList = [infile]
            )

MCTreeName="I3MCTree_postprop"
photonSeriesName = "aname"

kwargs = {}

tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
    GCDFile = options.gcdfile,
    DOMRadius=0.21590*icetray.I3Units.m, # 13" diameter 
    #IceModelLocation=options.ICEMODEL,
    IceModelLocation=Medium.MakePoneMediumProperties(),
    WavelengthAcceptance = dom_properties.GetCLSimQETable( factor=dom_properties.GetMaxAngularAcceptance()*1.05 ),
    MCPESeriesName = "",
    #MCPESeriesName = "MCPESeriesMap",
    
    PhotonSeriesName = photonSeriesName,
    MCTreeName = MCTreeName,
    RandomService = randomService,
    DOMEfficiency = 0.95,
    UseGPUs=not options.USECPU,
    UseCPUs=options.USECPU,
    #UseOnlyDeviceNumber=options.DEVICE,
    #UseCUDA=options.CUDA,
    EnableDoubleBuffering=options.double_buffering,
    UseI3PropagatorService=False,
    DOMOversizeFactor=options.OVERSIZE,
    UnWeightedPhotons=options.unweighted_photons,
    #DOMRadius = (17.0*2.54*0.01/2.0)*icetray.I3Units.m,
    CrossoverEnergyEM=None,
    PhotonHistoryEntries=0,
                #CrossoverEnergyHadron=float(options.CROSSENERGY),
    StopDetectedPhotons=True,
                #UseHoleIceParameterization=False, # Apply it when making hits!
    HoleIceParameterization=os.getenv('PONESRCDIR')+"/data/as.full",
    DoNotParallelize=False,
    UnshadowedFraction=1., #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
    )

tray.AddModule("I3Writer","writer",
#               SkipKeys=SkipKeys,
               Filename =  outfile,
               Streams = [icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation, icetray.I3Frame.DAQ],
              )


tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()

