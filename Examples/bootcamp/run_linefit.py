from icecube import icetray, phys_services, dataio
from icecube.icetray import I3Tray
from icecube.hdfwriter import I3HDFWriter
import icecube
icecube.icetray.load("linefit", False)

# import linefit

infile = '/project/6061446/cmiller/example_data/jp_linefit.i3.zst'
outfile = '/scratch/jpyanez/test_newlinefit.i3.zst'
hdf_filename = '/scratch/jpyanez/test.hdf'
all_i3_files = [infile, infile]

tray = I3Tray()

tray.AddModule("I3Reader", "i3_file_reader",
               FilenameList=all_i3_files)

## Add linefit here
tray.AddModule("I3LineFit","New_linefit_final_why_not",
               Name = "NEW_LINEFIT",
               InputRecoPulses = "PMTResponse_nonoise",
               LeadingEdge= "FLE",
               AmpWeightPower=0.0)

tray.Add(I3HDFWriter, 'MyWriter',
         Output=hdf_filename,
         Keys = ['JPLineFit'],
         SubEventStreams = [""])

# Add the I3File writer here
tray.AddModule("I3Writer","writer",
               Filename = outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
               )

tray.Execute()
tray.Finish()
