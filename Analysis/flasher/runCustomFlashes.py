from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
import genCustomFlashes
import genIsotropicFlashes
import genIsotropicFlashesSigma
from I3Tray import I3Units

from I3Tray import *
from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services

import math
import numpy
from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_flashesSC.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numframes", type="int", default=1,
                  dest="NUMFRAMES", help="The number of events per run")
parser.add_option("-c", "--CandleNumber", type="int", default=2,
                  dest="CANDLENUMBER", help="The number of events per run")
parser.add_option("-p", "--CandlePosition", type="float", default=107.66,
                  dest="CANDLEPOSITION", help="The number of events per run")

(options,args) = parser.parse_args()

if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = 12344,
    nstreams = 10000,
    streamnum = int(options.RUNNUMBER))

tray.context['I3RandomService'] = randomService

tray.AddModule("I3InfiniteSource","streams",
               Prefix='/home/users/akatil/P-ONE/GCD_files/STRAW_DOM_at_zero.i3.gz',  #STRAW_DOM_at_zero.i3.gz', #STRAW_sDOM_right.i3.gz', #STRAW_DirectDetection.i3.gz', #P_ONE_geometry_5_10_zeroPointOne_calib_updated.i3.gz', #STRAW_geometry.i3.gz
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               #Year=2012,
               #DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

tray.AddModule(genIsotropicFlashesSigma.CustomFlashes, "customFlasher",
		       FlasherPulseSeriesName = "CustomFlashes",
               PhotonsPerPulse = 1.7e6,
               CandleNumber = options.CANDLENUMBER,
               CandlePosition = options.CANDLEPOSITION)

tray.AddModule("I3Writer","writer",
    Filename = str(options.OUTFILE),
    Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

tray.Execute(options.NUMFRAMES+3)
