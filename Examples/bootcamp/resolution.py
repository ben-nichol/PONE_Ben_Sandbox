from icecube import icetray, phys_services, dataio, dataclasses
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


### Add a module that operates on I3 file

# Option 1 - just write a python function
def MyFunction(frame, RecoKey, ResultKey='RecoAngleError'):
    primary = frame['I3MCTree'].get_primaries()[0]
    reco_particle = frame[RecoKey]
    angle = phys_services.I3Calculator.angle(primary, reco_particle)
    print(angle)
    frame[ResultKey] = dataclasses.I3Double(angle)
    return True

tray.AddModule(MyFunction, 'calculate_resolution',
               RecoKey = 'JPLineFit',
               ResultKey = 'RecoAngleError')


tray.Add(I3HDFWriter, 'MyWriter',
         Output=hdf_filename,
         Keys = ['JPLineFit'],
         SubEventStreams = [""])

# Add the I3File writer here
tray.AddModule("I3Writer","writer",
               Filename = outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
               )

tray.Execute(40)
tray.Finish()
