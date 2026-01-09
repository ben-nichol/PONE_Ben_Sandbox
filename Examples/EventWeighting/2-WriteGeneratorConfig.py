from icecube import icetray, dataio, hdfwriter, LeptonInjector

tray = icetray.I3Tray()

# Read in a file that contains the S frame with LI config
tray.AddModule("I3Reader", "reader",
    Filename="muons.i3"
)

# Serialize the LI config from the S frame into a .lic file
tray.AddModule("InjectionConfigSerializer", "lic_writer",
    OutputPath="config.lic"
)

# Write Q frame data to hdf5 files
tray.AddSegment(hdfwriter.I3SimHDFWriter, "hdf",
    Output="output.h5",
    Keys=[
        "I3EventHeader",
        "EventProperties",
        "I3MCTree",
    ]
    )

tray.Execute()

