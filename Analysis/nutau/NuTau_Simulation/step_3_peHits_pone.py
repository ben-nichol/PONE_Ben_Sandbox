#!/usr/bin/env python

import sys
sys.path.insert(0,'/home/users/dhilu/P_ONE_dvirhilu/src')

from icecube import dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from experimentModelCode import FunctionClasses
from simAnalysis.SimAnalysis import passFrame
from os.path import expandvars
import numpy as np
import argparse
#from RemoveLatePhotons_V5 import RemoveLatePhotons

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-n', '--runNum',  dest = 'runNum', help = "number assigned to this specific run", default = 0 )
#parser.add_argument('-s', '--simType', dest = 'simType', help="which sim tool is used?")
parser.add_argument('-g', '--gcdType', dest = 'GCDType', help = "the type of GCD File used in the simulation")
parser.add_argument('-d', '--domType', dest = 'DOMType', help = "the type of DOM used in the simulation")
parser.add_argument('-H', '--hitThresh', help = "number of total hits required to not cut a frame")
parser.add_argument('-D', '--domThresh', help = "Number of DOMs with hits required to not cut a frame")
parser.add_argument('-f', '--filePath', help = "path of files will depend on whether code is run locally or not")
parser.add_argument('-i', '--infile', dest = 'infile', help= 'input file and path')
parser.add_argument('-o', '--outfile', dest = 'outfile', help= 'output file and path')
args = parser.parse_args()

if args.GCDType == 'testString':
    gcdPath = str(args.filePath) + 'I3Files/gcd/testStrings/HorizTestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = str(args.filePath) + 'I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz'
elif args.GCDType == 'denseGeo':
    gcdPath = str(args.filePath) + 'I3Files/gcd/denseGeo/denseGeo_n30_b50.0_a4.5_l7_linear_reset_offset_simple_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = str(args.filePath) + 'I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
elif args.GCDType == 'cube':
    gcdPath = str(args.filePath) + 'I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz'
elif args.GCDType == 'PONE':
    gcdPath = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")

#outname = 'nugen/nugenStep3/' + str(args.GCDType) + '/NuGen_step3_' + str(args.GCDType) + '_' + str(args.runNum) + '.i3.gz'
outname = 'genHitsTrial.i3.gz'
#inPath = str(args.filePath) + 'I3Files/nugen/nugenStep2/'+ str(args.GCDType) + '/NuGen_step2_' + str(args.GCDType) + '_' + str(args.runNum) + '.i3.gz'
inPath = args.infile
outPath = args.outfile

infile = dataio.I3File(inPath)
geofile = dataio.I3File(gcdPath)
outfile = dataio.I3File(outPath, 'w')

# get DOM characteristics
inFolder = '/home/users/akatil/P-ONE/git/PONE_NuTau/DOM/' + args.DOMType + '/'
filenameDomEff = 'DOMEfficiency.dat'
filenameAngAcc = 'AngularAcceptance.dat'

angAcc = FunctionClasses.Polynomial(np.loadtxt(inFolder + filenameAngAcc, ndmin = 1), -1, 1)
domEffData = np.loadtxt(inFolder + filenameDomEff, unpack = True)
domEff = FunctionClasses.FunctionFromTable(domEffData[0], domEffData[1])

# get geometry
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
calibration = cframe["I3Calibration"]
calMap = calibration.dom_cal

DOMRadius = 0.16510*I3Units.m
DOMOversizeFactor = 5.0
icemodel_efficiency_factor = 1.0
UnshadowedFraction = 1.0
HoleIceParameterization = expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm")
domAcceptance = clsim.GetIceCubeDOMAcceptance(
	                                domRadius = DOMRadius*DOMOversizeFactor,
	                                efficiency=icemodel_efficiency_factor*UnshadowedFraction)
domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=HoleIceParameterization)

def getSurvivalProbability(photon, omkey):
    if omkey in calMap:
        domcal = calMap[omkey]
        relativeDOMEff = domcal.relative_dom_eff * domcal.combined_spe_charge_distribution.compensation_factor
    else:
        relativeDOMEff = 1

    domGeo = geoMap[omkey]
    domDirection = domGeo.direction
    photonDirection = photon.dir
    dotProduct = photonDirection.x*domDirection.x + photonDirection.y*domDirection.y + photonDirection.z*domDirection.z
    # photon coming in, direction coming out. Sign on dot product should be flipped
    # directions are unit vectors already so cos_theta = -dotProduct (due to sign flip)
    probAngAcc = angAcc.getValue(-dotProduct)

    probDOMAcc = domEff.getValue(photon.wavelength)

    #print photon.weight, "Photon Weight"

    return probAngAcc*probDOMAcc*photon.weight*relativeDOMEff

def survived(photon,omkey):
    probability = getSurvivalProbability(photon, omkey)
    randomNumber = np.random.uniform()

    print "probability - ", probability

    if(probability > 1):
        #probability = ProbabilityCorrection(photon, omkey)
        #print "probability - ", probability
        raise ValueError("Probability = " + str(probability) + " > 1. Most likely photon weights are too high")

    if(probability > randomNumber):
        return True

    return False

def generateMCPEList(photons, modkey):
    omkey = OMKey(modkey.string, modkey.om, 0)
    mcpeList = simclasses.I3MCPESeries()
    photonList = []
    for photon in photons:
        if survived(photon,omkey):
            mcpe = simclasses.I3MCPE()
            #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
            mcpe.npe = 1
            mcpe.time = photon.time #TODO: change to corrected time
            mcpeList.append(mcpe)
            photonList.append(photon)

    return mcpeList, photonList

# TODO: def getCorrectedTime(photon, omkey):

while( infile.more() ):
    frame = infile.pop_daq()
    photonDOMMap = frame["I3Photons"]
    mcpeMap = simclasses.I3MCPESeriesMap()
    succPhotonMap = simclasses.I3CompressedPhotonSeriesMap()
    for modkey in photonDOMMap.keys():
        mcpeList, photonList = generateMCPEList(photonDOMMap[modkey], modkey)
        omkey = OMKey(modkey.string, modkey.om, 0)
        if len(mcpeList) > 0:
            mcpeMap[omkey] = mcpeList
            succPhotonMap[modkey] = photonList

    if len(mcpeMap) > 0:
        frame["MCPESeriesMap"] = mcpeMap

    # only add frame to file if a hit was generated
    #if passFrame(frame, mcpeMap.keys(), int(args.hitThresh), int(args.domThresh), ):
        frame["SuccPhotonMap"] = succPhotonMap
        #frame.Delete("I3Photons")
        outfile.push(frame)

outfile.close()
