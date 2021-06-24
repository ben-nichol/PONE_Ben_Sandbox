from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
import numpy as np
import random, re

'''
Functions for MCPEConverter.py
MCPEConverter.py makes hits from I3Photons
'''
def getSurvivalProbability(geoMap, photon, omkey, angularAcceptance, domAcceptance, relativeEfficieny):

    domGeo = geoMap[omkey]
    domDirection = domGeo.direction
    photonDirection = photon.dir
    dotProduct = photonDirection.x*domDirection.x + photonDirection.y*domDirection.y + photonDirection.z*domDirection.z
    # photon coming in, direction coming out. Sign on dot product should be flipped
    # directions are unit vectors already so cos_theta = -dotProduct (due to sign flip)

    probDOMAcc = domAcceptance.GetValue(photon.wavelength)

    probAngAcc = angularAcceptance

    return probAngAcc*probDOMAcc*photon.weight*relativeEfficieny


def survived(geoMap, photon,omkey, angularAcceptance, domAcceptance, relativeEfficieny):
    probability = getSurvivalProbability(geoMap, photon, omkey, angularAcceptance, domAcceptance, relativeEfficieny)
    randomNumber = np.random.uniform()

    if(probability > 1):

        raise ValueError("Probability = " + str(probability) + " > 1. Most likely photon weights are too high")

    if(probability > randomNumber):
        return True

    return False

def generateMCPEList(geoMap, photons, modkey, angularAcceptance, domAcceptance, relativeEfficieny):
    omkey = OMKey(modkey.string, modkey.om, 0)
    mcpeList = simclasses.I3MCPESeries()
    photonList = []
    for photon in photons:
        if survived(geoMap,photon,omkey,angularAcceptance,domAcceptance,relativeEfficieny):
            mcpe = simclasses.I3MCPE()
            #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
            mcpe.npe = 1
            mcpe.time = photon.time #+ np.random.normal(0, 1.5) #smearing the time with random values from a gaussian with width 1.5ns
            mcpeList.append(mcpe)
            photonList.append(photon)

    return mcpeList, photonList

'''
Functions for NoiseGenerator.py
NoiseGenrator.py injects noise hits from the STRAW data
'''

def selectStartTime(timestamps):
    randStartTime = random.uniform(0, max(timestamps))
    noiseWindow = timestamps[(timestamps >= randStartTime) & (timestamps <= randStartTime+72e2)]
    return randStartTime, noiseWindow

def addHits(timestamps):
    print "adding Hits"
    mDOM_noise = ([])
    startTimeVals = ([])
    for segment in range(0, 24):
        startTime, noiseWindow = selectStartTime(timestamps)
        movingVals = (noiseWindow - startTime)-(72e2/2)
        startTimeVals = np.append(startTimeVals, startTime)
        mDOM_noise = np.append(mDOM_noise, movingVals)
    #print "mDOM_noise", mDOM_noise
    return mDOM_noise

def generateNoiseMCPEList(modkey, medVal, timestamps):
    omkey = OMKey(modkey.string, modkey.om, 0)
    noiseHits = addHits(timestamps)
    noiseList = simclasses.I3MCPESeries()

    for noiseHit in noiseHits:
        mcpeNoise = simclasses.I3MCPE()
        #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
        mcpeNoise.npe = 1
        mcpeNoise.time = noiseHit + medVal
        #print "mcpeNoise.time", mcpeNoise.time
        noiseList.append(mcpeNoise)

    #print 'NoiseList', noiseList
    return noiseList

'''
Functions for RecoPulseGenerator.py
Creating RecoPulses from MCPEs
'''

def generateRecoPulses(omkey, hits):
    recoPulses = dataclasses.I3RecoPulseSeries()
    for mcpe in hits:
        recoPulse = dataclasse.I3RecoPulse()
        recoPulse.time = mcpe.time + np.random.normal(0, 1.5) #smearing the time with random values from a gaussian with width 1.5ns
        recoPulse.charge = mcpe.npe
        #recoPulse.width = (Leave it blank for now!!!)

        recoPulses.append(recoPulse)

    return recoPulses
