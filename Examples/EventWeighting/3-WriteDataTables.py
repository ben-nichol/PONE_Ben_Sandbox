from icecube import icetray, dataio, hdfwriter, dataclasses, LeptonInjector
from icecube.icetray import I3Tray

tray = I3Tray()

tray.AddModule("I3Reader", "reader",
    Filename="muons.i3"
)

tray.AddSegment(hdfwriter.I3SimHDFWriter, "hdf",
    Output="output.h5",
    Keys=[
        "I3EventHeader",
        "EventProperties",
        "I3MCTree",
    ]
    )

tray.Execute()

