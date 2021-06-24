from I3Tray import *
from icecube import icetray, dataclasses, phys_services, sim_services, dataio,  earthmodel_service, neutrino_generator, tableio, hdfwriter
from icecube.simprod import segments
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
from icecube.simclasses import I3MMCTrack
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "A scripts to run the neutrino generation simulation step using Neutrino Generator")

parser.add_argument('-emin', '--energyMin', default = 5.0, help = "the minimum energy")
parser.add_argument('-emax', '--energyMax', default = 8.0, help = "the maximum energy")
parser.add_argument('-n', '--numEvents', help = "number of events produced by the simulation")
parser.add_argument('-o', '--outfile', help="name and path of output file")
parser.add_argument('-r', '--runNum', help="run Number")

args = parser.parse_args()
emin = float(args.energyMin)
emax = float(args.energyMax)
print(emin, emax)
numEvents = int(args.numEvents)
runNum = int(args.runNum)

cylinder = [float(500), float(1000), float(0), float(0), float(0)]
zenithMin = (90 - float(45)) * I3Units.deg
zenithMax = (90 + float(45)) * I3Units.deg
azimuthMin = 0 * I3Units.deg
azimuthMax = 180 * I3Units.deg
#gcd = "/home/users/jguthrie/oscnext/GCD_files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz"
gcd = "/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz"

tray = I3Tray()

args_ratios = "1:1"
ratio = [float(ratio) for ratio in args_ratios.split(":")]

#Random
#randomService  = phys_services.I3GSLRandomService(123456)
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
               NuTypes = ['NuTau','NuTauBar'],
               PrimaryTypeRatio = ratio,
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
                TablesDir = "/home/users/dhilu/P_ONE_dvirhilu/CrossSectionModels",
                CrossSectionModel = "csms_differential_v1.0"
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

from icecube import PROPOSAL
from os.path import expandvars

#propagators = sim_services.I3ParticleTypePropagatorServiceMap()

#mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/mediadef')

#TauMinusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
		#mediadef=mediadef,
		#cylinderRadius=500,
		#cylinderHeight=1000,
		#type=I3Particle.ParticleType.TauMinus)

#TauMinusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
		#config_file='/home/users/akatil/software/V06-01-02/build/PROPOSAL/resources/config_icesim.json',
		#config_file='/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/V00-00-02/PROPOSAL/resources/config_icesim.json',
		#final_stochastic_loss=I3Particle.ParticleType.TauMinus,
		#distance_to_propagate = 1e+20)

#TauPlusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
		#config_file='/home/users/akatil/software/V06-01-02/build/PROPOSAL/resources/config_icesim.json',
		#config_file='/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/V00-00-02/PROPOSAL/resources/config_icesim.json',
		#final_stochastic_loss=I3Particle.ParticleType.TauPlus,
		#distance_to_propagate = 1e+20)

#propagators[I3Particle.ParticleType.TauMinus] = TauMinusPropagator
#propagators[I3Particle.ParticleType.TauPlus] = TauPlusPropagator

#tray.Add('I3PropagatorModule', 'Tau_propagator',
		#PropagatorServices=propagators,
		#RandomService=randomService,
		#InputMCTreeName="I3MCTree_NuGen",
        #OutputMCTreeName="I3MCTree")

tray.Add(segments.PropagateMuons, 'ParticlePropagators',
         RandomService=randomService,
         SaveState=True,
         InputMCTreeName="I3MCTree_NuGen",
         OutputMCTreeName="I3MCTree")

#converts q frames?
tray.Add("I3NullSplitter",
       SubEventStreamName = "fullevent")

#creates HDF5 file
#tray.Add(tableio.I3TableWriter,'writer',
        #ableservice=hdfwriter.I3HDFTableService(args.outfile+'.hdf5'),
	 	#keys=["I3MCTree"],
        #SubEventStreams=["fullevent"])

#tray.Add("I3Writer", filename ='NuTau_test_energyRanges.i3',
#        streams = [icetray.I3Frame.DAQ],)

tray.Add("I3Writer", filename = args.outfile+'.i3.gz',
        streams = [icetray.I3Frame.DAQ],)

tray.Execute()
