from I3Tray import *
from icecube import icetray, dataclasses, phys_services, sim_services, dataio,  earthmodel_service, neutrino_generator, tableio, hdfwriter
from icecube.simprod import segments
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
from icecube.simclasses import I3MMCTrack
from icecube import PROPOSAL
import os
from os.path import expandvars
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "A scripts to run the neutrino generation simulation step using Neutrino Generator")

parser.add_argument('-emin', '--energyMin', default = 5.0,                                            help="the minimum energy")
parser.add_argument('-emax', '--energyMax', default = 8.0,                                            help="the maximum energy")
parser.add_argument('-n',    '--numEvents', default = 1000,                                           help="number of events produced by the simulation")
parser.add_argument('-o',    '--outfile',   default = "output.i3",                                    help="name and path of output file")
parser.add_argument('-r',    '--runNum',    default = 0,                                              help="run Number")
parser.add_argument("-a",    "--ratios",    default="1.0:1.0:1.0:1.0:1.0:1.0",                                help="ratio of input neutrino")
parser.add_argument("-t",    "--types",     default="NuE:NuEBar:NuTau:NuTauBar:NuMu:NuMuBar",                      help="type of input neutrino")
parser.add_argument("-g",    "--gcd",       default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz", help="gdc file")
parser.add_argument("-c",    "--crossdir",  default=os.getenv('PONESRCDIR')+"/CrossSectionModels",    help='path to cross section models')
parser.add_argument("-m",    "--crossmodel",default='csms_differential_v1.0',                         help='cross section model')

args = parser.parse_args()

typeString = args.types
ratioString = args.ratios

typevec = typeString.split(":")
ratiostvec = ratioString.split(":")
ratiovec = []
for ratio in ratiostvec:
    ratiovec.append(float(ratio))

emin = float(args.energyMin)
emax = float(args.energyMax)
numEvents = int(args.numEvents)
runNum = int(args.runNum)
print(emin, emax, ratiovec, typevec, numEvents, runNum)

cylinder = [float(300), float(1300), float(0), float(0), float(0)]
zenithMin = 0 * I3Units.deg
zenithMax = 180 * I3Units.deg
azimuthMin = 0 * I3Units.deg
azimuthMax = 360 * I3Units.deg

gcd = args.gcd

tray = I3Tray()

#Random
randomService = phys_services.I3SPRNGRandomService(
        seed = 1234657,
		nstreams = 10000,
        streamnum = runNum)
tray.context['I3RandomService'] = randomService
tray.Add("I3InfiniteSource", prefix = gcd)

tray.Add("I3MCEventHeaderGenerator",
	       EventID=1,
	       IncrementEventID=True)

tray.Add("I3EarthModelServiceFactory", "EarthModelService",
                EarthModels = ["PREM_mmc"],
                MaterialModels = ["Standard"],
                IceCapType = "IceSheet",
                DetectorDepth = 2600*I3Units.m,
                PathToDataFileDir = "")

tray.Add("I3NuGSteeringFactory", "steering",
                EarthModelName = "EarthModelService",
                NEvents = numEvents,
                SimMode = "Detector",
                VTXGenMode = "NuGen",
                InjectionMode = "surface",
                CylinderParams = cylinder,
                DoMuonRangeExtension = False,
                UseSimpleScatterForm = True,
                MCTreeName = "I3MCTree_NuGen"
                )

tray.Add("I3NuGDiffuseSource","diffusesource",
               RandomService = randomService,
               SteeringName = "steering",
		#NuFlavor = 'NuTau',
               NuTypes = typevec,#['NuTau','NuTauBar'],
               PrimaryTypeRatio = ratiovec,
               GammaIndex = 2.19,
               EnergyMinLog = emin,
               EnergyMaxLog = emax,
               ZenithMin = zenithMin,
               ZenithMax = zenithMax,
               AzimuthMin = azimuthMin,
               AzimuthMax = azimuthMax,
               ZenithWeightParam = 1.0,
               AngleSamplingMode = "COS"
              )

tray.Add("I3NuGInteractionInfoDifferentialFactory", "interaction",
                RandomService = randomService,
                SteeringName = "steering",
                TablesDir = args.crossdir,
                CrossSectionModel = args.crossmodel
              )

tray.Add("I3NeutrinoGenerator","generator",
                RandomService = randomService,
                SteeringName = "steering",
                InjectorName = "diffusesource",
                InteractionInfoName = "interaction",
                #PropagationWeightMode = "NCGRWEIGHTED",
                InteractionCCFactor = 1.0,
                InteractionNCFactor = 1.0,
                #InteractionGRFactor = 1.0
              )

tray.Add(segments.PropagateMuons, 'ParticlePropagators',
         RandomService=randomService,
         SaveState=True,
         InputMCTreeName="I3MCTree_NuGen",
         OutputMCTreeName="I3MCTree")

#converts q frames?
tray.Add("I3NullSplitter",
       SubEventStreamName = "fullevent")

tray.Add("I3Writer", filename = args.outfile,
        streams = [icetray.I3Frame.DAQ],)

tray.Execute()
