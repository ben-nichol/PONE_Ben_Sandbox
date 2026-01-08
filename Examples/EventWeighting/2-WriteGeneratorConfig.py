from icecube import icetray, dataio, LeptonInjector

tray = icetray.I3Tray()

# Read in a file that contains the S frame with LI config
tray.AddModule("I3Reader", "reader",
    Filename="muons.i3"
)

# Serialize the LI config from the S frame into a .lic file
tray.AddModule("InjectionConfigSerializer", "lic_writer",
    OutputPath="config.lic"
)

tray.Execute()

