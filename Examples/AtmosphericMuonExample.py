# Creates 10 atmospheric muons above the detector and propagates them.

from I3Tray import I3Tray
from icecube import icetray, phys_services
from icecube.icetray import I3Units
from segments import PropagateMuons
from Generators import AtmosphericMuons

rand = phys_services.I3SPRNGRandomService(
    seed = 0, 
    nstreams = 100000000, 
    streamnum = 0)

tray = I3Tray()
tray.Add("I3InfiniteSource")
tray.AddModule(AtmosphericMuons.MuonGenerator,
    RandomService=rand)
tray.Add(PropagateMuons,
    RandomService=rand,
    SaveState=True,
    InputMCTreeName="MuonGeneratorI3MCTree",
    OutputMCTreeName="muon_track")
tray.AddModule("I3Writer",
    Filename="AtmosphericMuonExample.i3",
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Simulation])
tray.AddModule("TrashCan", "adios")
tray.Execute(10)
tray.Finish()
