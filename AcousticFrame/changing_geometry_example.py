from icecube import icetray, dataclasses, phys_services, clsim
from AcousticFrame import geometry, frames
import os, numpy as nn

tray = icetray.I3Tray()
tray.AddModule("I3InfiniteSource", Stream=icetray.I3Frame.DAQ)

def createEvents(frame):
    energy    = 1e3 * icetray.I3Units.GeV
    position  = dataclasses.I3Position(0,10,500)
    direction = dataclasses.I3Direction(0.,0.,-1.)
    primary = dataclasses.I3Particle()
    primary.type = dataclasses.I3Particle.ParticleType.NuE
    primary.location_type = dataclasses.I3Particle.LocationType.Anywhere
    daughter = dataclasses.I3Particle()
    daughter.type = dataclasses.I3Particle.ParticleType.EMinus
    daughter.location_type = dataclasses.I3Particle.LocationType.InIce
    for particle in [primary, daughter]:
        particle.energy = energy
        particle.pos = position
        particle.dir = direction
        particle.time = 0.
    mctree = dataclasses.I3MCTree()
    mctree.add_primary(primary)
    mctree.append_child(primary,daughter)
    frame["I3MCTree"] = mctree
tray.Add(createEvents, Streams=[icetray.I3Frame.DAQ])

cm = icetray.I3Units.centimeter
s  = icetray.I3Units.s
for delay, speed, angle in [[ 0, 20*cm/s, 0],
                            [ 5, 10*cm/s, 0],
                            [10,  0*cm/s, 0],
                            [15, 10*cm/s, nn.pi],
                            [20, 20*cm/s, nn.pi]]:
    tray.AddModule(frames.injectAcousticFrame,
        delay=delay, current_speed=speed, current_angle=angle)
tray.AddModule(frames.injectDetectorStatusFrame,
    insertbefore=frames.Acoustic)
tray.AddModule(frames.injectCalibrationFrame,
    insertbefore=frames.Acoustic)

tray.AddModule(frames.createGeometryFromAcoustic)
tray.AddModule("I3Writer", filename='changing_geometry_example.i3')

tray.AddModule(frames.interpolateGeometryFromAcoustic)
tray.AddModule("I3Writer", filename='changing_geometry_example_interpolated.i3')

randomService = phys_services.I3SPRNGRandomService(seed=0, nstreams=1, streamnum=0)
class WrapCLSim(frames.CreateSeparateGCD):
    def InsertModule(self, wrapper_tray, gcdfile):
        # This is just testing CLSim within the wrapper. Do not use these CLSim settings for P-ONE simulations.
        wrapper_tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
            GCDFile=gcdfile,
            PhotonSeriesName='PhotonSeries',
            MCTreeName="I3MCTree",
            RandomService=randomService,
            UseGPUs=True,
            UseCPUs=False,
            UseOnlyDeviceNumber=None,
            UseCUDA=False,
            EnableDoubleBuffering=False,
            UseI3PropagatorService=False,
            IceModelLocation=os.path.expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_lea"),
            UnWeightedPhotons=False)
tray.AddModule(WrapCLSim)
def checkPhotons(frame):
    n = sum([len(p) for k,p in frame['PhotonSeries']])
    print("%i photons detected" % n)
tray.Add(checkPhotons, Streams=[icetray.I3Frame.DAQ])
tray.AddModule("I3Writer", filename='changing_geometry_example_clsim.i3')

tray.Execute(21)
tray.Finish()

