from icecube import icetray, phys_services, dataio
from icecube.icetray import I3Tray
from icecube.hdfwriter import I3HDFWriter


infile = '/project/6061446/cmiller/example_data/jp_linefit.i3.zst'
hdf_filename = '/scratch/jpyanez/test.hdf'
all_i3_files = [infile, infile]

tray = I3Tray()

tray.AddModule("I3Reader", "i3_file_reader",
               FilenameList=all_i3_files)

tray.Add(I3HDFWriter, 'MyWriter',
         Output=hdf_filename,
         Keys = ['JPLineFit'],
         SubEventStreams = [""])

tray.Execute()
tray.Finish()
