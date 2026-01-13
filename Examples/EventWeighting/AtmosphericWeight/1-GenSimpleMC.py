from icecube import icetray, dataclasses, phys_services, dataio
from icecube.icetray import I3Tray, I3Units
from icecube import LeptonInjector

tray = I3Tray()

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

injector_list = []
injector_list.append(
    LeptonInjector.injector(
        NEvents=10,
        FinalType1=dataclasses.I3Particle.ParticleType.MuMinus,
        FinalType2=dataclasses.I3Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile="/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/dsdxdy_nu_CC_iso.fits",
        TotalCrossSectionFile="/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/sigma_nu_CC_iso.fits",
        Ranged=True,
    )
    )
    # Create the multileptoninjector object with your list of injectors
tray.AddModule(
    "MultiLeptonInjector",
    EarthModel="Earth",
    Generators=injector_list,
    MinimumEnergy=(10.0 ** (3)) * I3Units.GeV,
    MaximumEnergy=(10.0 ** (4)) * I3Units.GeV,
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


tray.Add("I3Writer",filename="muons.i3")
tray.Execute()
