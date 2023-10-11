

import argparse
from os.path import expandvars
import os, sys, random
import random
from icecube import icetray, dataclasses, dataio, simclasses
from icecube import phys_services, sim_services
from icecube import clsim
from icecube.icetray import I3Tray
import WaterOpticalModel.MakePoneMediumPropertiesConservative as Medium
#import WaterOpticalModel.MakePoneMediumPropertiesSpeculativeExtendedRange as Medium
from Utilities.DOMUtility import DOMProperties

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument("-i", "--infile",default="dataio/muonprop.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-o", "--outfile",default="dataio/photonemit.i3", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-g", "--gcdfile",default=os.getenv('PONESRCDIR')+"/GCD/PONE_5String.i3.gz", help="Readin GCD file")
args = parser.parse_args()
count = 0
CPU=False

# load DOM properties
dom_properties = DOMProperties()

photon_series = "I3Photons"
#print 'CUDA devices: ', options.DEVICE
tray = I3Tray()

# Now fire up the random number generator with that seed
#from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
                seed = 0,
                nstreams = int(4e7),
                streamnum = 0)

tray.context['I3RandomService'] = randomService

outfile = args.outfile

infile = args.infile

#gcd_file = dataio.I3File(args.gcdfile)
print(args.gcdfile)

tray.AddModule('I3Reader', 'reader',
            FilenameList = [infile,args.gcdfile]
            )

print("GetMaxTotalAcceptance()   =", dom_properties.GetMaxTotalAcceptance())
print("GetMaxAngularAcceptance() =", dom_properties.GetMaxAngularAcceptance())
print("GetMaxPMTQE()             =", dom_properties.GetMaxPMTQE())

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseCPUs=True,
#                UseGPUs=True,
                #UseOnlyDeviceNumber=[1],
                #OpenCLDeviceList=[0],
                MCTreeName="I3MCTree",
                UseI3PropagatorService=False,
                #OutputMCTreeName="I3MCTree_clsim",
                #FlasherInfoVectName="I3FlasherInfo",
                #FlasherPulseSeriesName="PhotonBomb",
                #MMCTrackListName="MMCTrackList",
                PhotonSeriesName=photon_series,
                MCPESeriesName='',
                #ParallelEvents=1000,
                RandomService=randomService,
                IceModelLocation=Medium.MakePoneMediumProperties(),
                UnWeightedPhotons=False,
                #UnWeightedPhotonsScalingFactor = None,
		DOMRadius = (17.0*2.54*0.01/2.0)*icetray.I3Units.m,
                UseGeant4=False,
                CrossoverEnergyEM=None,
                PhotonHistoryEntries=0,
                #CrossoverEnergyHadron=float(options.CROSSENERGY),
                StopDetectedPhotons=True,
                #UseHoleIceParameterization=False, # Apply it when making hits!
                #HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.flasher_p1_0.30_p2_-1"),
                DoNotParallelize=False,
		WavelengthAcceptance = dom_properties.GetCLSimQETable( factor=dom_properties.GetMaxAngularAcceptance()*1.05 ),
                DOMOversizeFactor=1.0, #(17./13.),
                UnshadowedFraction=1., #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile= args.gcdfile, #gcd_file,
                )

#icetray.logging.I3Logger.global_logger.set_level_for_unit('clsim', icetray.logging.I3LogLevel.LOG_ERROR)
#icetray.logging.I3Logger.global_logger.set_level_for_unit('I3CLSimStepToPhotonConverterOpenCL', icetray.logging.I3LogLevel.LOG_WARN)

#tray.AddModule(PrintMessage,"print",message = "CLSiM Check")

tray.AddModule("I3Writer","writer",
#               SkipKeys=SkipKeys,
               Filename =  outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

#tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
