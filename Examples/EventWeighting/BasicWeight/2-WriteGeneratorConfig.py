from icecube import icetray, dataio, hdfwriter, LeptonInjector
import sys

if len(sys.argv) < 3:
    print("Usage: {} <input1.i3> [<input2.i3> ...] <output>".format(sys.argv[0]))
    sys.exit(1)

# All arguments except the last are input files
input_files = sys.argv[1:-1]
# The last argument is the output base name (for .h5 and .lic)
output = sys.argv[-1]

# Form the .lic file name by removing .h5 or .hdf5 and appending .lic
output_lic = output + '.lic'
output_h5 = output + '.h5'


tray = icetray.I3Tray()

# Read in all input files
tray.AddModule("I3Reader", "reader",
    FilenameList=input_files
)

# Serialize the LI config from the S frame into a .lic file
tray.AddModule("InjectionConfigSerializer", "lic_writer",
    OutputPath=output_lic
)

# Write Q frame data to hdf5 files
tray.AddSegment(hdfwriter.I3SimHDFWriter, "hdf",
    Output=output_h5,
    Keys=[
        "I3EventHeader",
        "EventProperties",
        "I3MCTree",
    ]
)

tray.Execute()

