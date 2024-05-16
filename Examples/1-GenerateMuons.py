from icecube import icetray, dataclasses, phys_services, dataio
from icecube.icetray import I3Tray, I3Units
from icecube import LeptonInjector
import os
import argparse

parser = argparse.ArgumentParser(
    description="A script to run the neutrino generation simulation step using Neutrino Generator"
)

parser.add_argument("-emin", "--energyMin", default=1.0, help="the minimum energy")
parser.add_argument("-emax", "--energyMax", default=3.0, help="the maximum energy")
parser.add_argument(
    "-n", "--numEvents", default=100, help="number of events produced by the simulation"
)
parser.add_argument(
    "-o", "--outfile", default="dataio/muons.i3", help="name and path of output file"
)
parser.add_argument(
    "-c",
    "--crossdir",
    default=os.getenv("PONESRCDIR") + "/CrossSectionModels/csms_differential_v1.0",
    help="path to cross section models",
)

args = parser.parse_args()


emin = float(args.energyMin)
emax = float(args.energyMax)
numEvents = int(args.numEvents)
print(emin, emax, numEvents)

tray = I3Tray()

# Random
randomService = phys_services.I3GSLRandomService(seed=12345)
tray.context["I3RandomService"] = randomService
tray.AddModule("I3InfiniteSource", "TheSource", Stream=icetray.I3Frame.DAQ)
tray.Add(
    "I3EarthModelServiceFactory",
    "Earth",
    EarthModels=["PREM_pone"],
    MaterialModels=["Standard"],
    IceCapType="IceSheet",
    DetectorDepth=(2600 - 500) * I3Units.m,
    PathToDataFileDir="",
)

# we create a list of injector objects
#   each of these injectors can have a unique cross sections used
injector_list = []
injector_list.append(
    LeptonInjector.injector(
        NEvents=args.numEvents,
        FinalType1=dataclasses.I3Particle.ParticleType.MuMinus,
        FinalType2=dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile=args.crossdir + "/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile=args.crossdir + "/sigma_nu_CC_iso.fits",
        Ranged=True,
    )
)
injector_list.append(
    LeptonInjector.injector(
        NEvents=args.numEvents,
        FinalType1=dataclasses.I3Particle.ParticleType.MuPlus,
        FinalType2=dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile=args.crossdir + "/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile=args.crossdir + "/sigma_nu_CC_iso.fits",
        Ranged=True,
    )
)

# Create the multileptoninjector object with your list of injectors
tray.AddModule(
    "MultiLeptonInjector",
    EarthModel="Earth",
    Generators=injector_list,
    MinimumEnergy=(10.0 ** (args.energyMin)) * I3Units.GeV,
    MaximumEnergy=(10.0 ** (args.energyMax)) * I3Units.GeV,
    MinimumZenith=0.0 * I3Units.deg,
    MaximumZenith=180.0 * I3Units.deg,
    PowerLawIndex=1.0,
    InjectionRadius=600 * I3Units.meter,
    EndcapLength=700 * I3Units.meter,
    CylinderRadius=700 * I3Units.meter,
    CylinderHeight=1000 * I3Units.meter,
    MinimumAzimuth=0.0 * I3Units.deg,
    MaximumAzimuth=360.0 * I3Units.deg,
    RandomService="I3RandomService",
)


tray.Add(
    "I3Writer",
    filename=args.outfile,
    streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ],
)

tray.Execute()
