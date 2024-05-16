from icecube import icetray, phys_services
from icecube.icetray import I3Tray
from segments import PropagateMuons
import os
import argparse

parser = argparse.ArgumentParser(description="Propagate muons through water")

parser.add_argument(
    "-i", "--infile", default="dataio/muons.i3", help="name and path of input file"
)
parser.add_argument(
    "-o", "--outfile", default="dataio/prop.i3", help="name and path of output file"
)

args = parser.parse_args()


tray = I3Tray()

# Random
randomService = phys_services.I3GSLRandomService(seed=4321)
tray.context["I3RandomService"] = randomService

tray.AddModule("I3Reader", "reader", FilenameList=[args.infile])

tray.Add(
    PropagateMuons,
    "ParticlePropagators",
    RandomService=randomService,
    SaveState=True,
    InputMCTreeName="I3MCTree",
    OutputMCTreeName="I3MCTree_postprop",
    PROPOSAL_config_file=os.getenv("PONESRCDIR") + "/configs/PROPOSAL_config.json",
)


tray.Add(
    "I3Writer",
    filename=args.outfile,
    streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ],
)

tray.Execute()
