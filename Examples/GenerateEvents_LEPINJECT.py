from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube import LeptonInjector 
from icecube.icetray import I3Units
from icecube import PROPOSAL
from segments import PropagateMuons
import os
from os.path import expandvars
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "A scripts to run the neutrino generation simulation step using Neutrino Generator")

parser.add_argument('-emin', '--energyMin', default = 3.0,                                            help="the minimum energy")
parser.add_argument('-emax', '--energyMax', default = 7.0,                                            help="the maximum energy")
parser.add_argument('-n',    '--numEvents', default = 1000,                                           help="number of events produced by the simulation")
parser.add_argument('-o',    '--outfile',   default = "output.i3",                                    help="name and path of output file")
parser.add_argument('-r',    '--runNum',    default = 0,                                              help="run Number")
parser.add_argument("-a",    "--ratios",    default="1.0:1.0:1.0:1.0:1.0:1.0",                                help="ratio of input neutrino")
parser.add_argument("-t",    "--types",     default="NuE:NuEBar:NuTau:NuTauBar:NuMu:NuMuBar",                      help="type of input neutrino")
parser.add_argument("-g",    "--gcd",       default=os.getenv('PONESRCDIR')+"/GCD/PONE_10String.i3.gz", help="gdc file")
parser.add_argument("-c",    "--crossdir",  default=os.getenv('PONESRCDIR')+"/CrossSectionModels/csms_differential_v1.0",    help='path to cross section models')
parser.add_argument("-x",    "--config",    default="",help="")

args = parser.parse_args()

typeString = args.types
ratioString = args.ratios

typevec = typeString.split(":")
ratiostvec = ratioString.split(":")
ratiovec = []
for ratio in ratiostvec:
    ratiovec.append(float(ratio))

emin = float(args.energyMin)
emax = float(args.energyMax)
numEvents = int(args.numEvents)
runNum = int(args.runNum)
print(emin, emax, ratiovec, typevec, numEvents, runNum)

cylinder = [float(300), float(1300), float(0), float(0), float(0)]
zenithMin = 0 * I3Units.deg
zenithMax = 180 * I3Units.deg
azimuthMin = 0 * I3Units.deg
azimuthMax = 360 * I3Units.deg

gcd = args.gcd

tray = I3Tray()

#Random
randomService = phys_services.I3GSLRandomService(seed=int(args.runNum))
tray.context["I3RandomService"] = randomService
tray.AddService("I3EarthModelServiceFactory", "Earth")
tray.AddModule("I3InfiniteSource", "TheSource", Stream=icetray.I3Frame.DAQ)

tray.Add("I3EarthModelServiceFactory", "Earth",
                EarthModels = ["PREM_mmc"],
                MaterialModels = ["Standard"],
                IceCapType = "IceSheet",
                DetectorDepth = 2600*I3Units.m,
                PathToDataFileDir = "")

# we create a list of injector objects
#   each of these injectors can have a unique cross sections used 
injector_list = []
injector_list.append(
    LeptonInjector.injector(
        NEvents         = args.numEvents,
        FinalType1      = dataclasses.I3Particle.ParticleType.MuMinus,
        FinalType2      = dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile  = args.crossdir + "/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile               = args.crossdir + "/sigma_nu_CC_iso.fits",
        Ranged = True)
    )
injector_list.append(
    LeptonInjector.injector(
        NEvents         = args.numEvents,
        FinalType1      = dataclasses.I3Particle.ParticleType.EMinus,
        FinalType2      = dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile  = args.crossdir + "/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile               = args.crossdir + "/sigma_nu_CC_iso.fits",
        Ranged = False)
    )
injector_list.append(
    LeptonInjector.injector(
        NEvents         = args.numEvents,
        FinalType1      = dataclasses.I3Particle.ParticleType.TauMinus,
        FinalType2      = dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile  = args.crossdir + "/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile               = args.crossdir + "/sigma_nu_CC_iso.fits",
        Ranged = False)
    )

# Create the multileptoninjector object with your list of injectors 
tray.AddModule("MultiLeptonInjector",
    EarthModel      = "Earth",
    Generators      = injector_list,
    MinimumEnergy   = (10.0**(args.energyMin)) * I3Units.GeV,
    MaximumEnergy   = (10.0**(args.energyMax)) * I3Units.GeV,
    MinimumZenith   = 90 * I3Units.deg,
    MaximumZenith   = 180 * I3Units.deg, 
    PowerLawIndex   = 2.,
    InjectionRadius = 500 * I3Units.meter,
    EndcapLength    = 1100 * I3Units.meter,
    CylinderRadius  = 100 * I3Units.meter,
    CylinderHeight  = 1000 * I3Units.meter,
    MinimumAzimuth  = 0. * I3Units.deg,
    MaximumAzimuth  = 360. * I3Units.deg,
    RandomService   = "I3RandomService")

tray.AddModule("InjectionConfigSerializer", OutputPath=args.config)

tray.Add(PropagateMuons, 'ParticlePropagators',
         RandomService=randomService,
         SaveState=True,
         InputMCTreeName="I3MCTree",
         OutputMCTreeName="I3MCTree_postprop",
         PROPOSAL_config_file=os.getenv('PONESRCDIR')+"/configs/PROPOSAL_config.json")

event_id = 1
def get_header(frame):
    global event_id 
    header          = dataclasses.I3EventHeader()
    header.event_id = event_id
    header.run_id   = int(args.runNum)
    frame["I3EventHeader"] = header

    event_id += 1

tray.AddModule(get_header, streams = [icetray.I3Frame.DAQ])

#converts q frames?
tray.Add("I3NullSplitter",
       SubEventStreamName = "fullevent")

tray.Add("I3Writer", filename = args.outfile+"_"+str(args.runNum)+".i3.gz",
        streams = [icetray.I3Frame.DAQ],)

tray.Execute()
