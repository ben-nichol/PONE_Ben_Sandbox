'''
-- SALT --

Module for generating K40 noise in P-ONE
'''

import numpy as np
import pickle as pickle
import os
from NoiseGenHelpers import makePhoton, unpickle, myQuaternion
from NoiseGenHelpers import K40

from icecube import icetray, dataio, dataclasses, simclasses
from icecube.icetray import OMKey, I3Units, I3Frame
from icecube.dataclasses import ModuleKey
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, DOMProperties

class K40Noise(icetray.I3ConditionalModule):

    def __init__(self, context):

        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("SEED","Random number generator seed",None)
        self.AddParameter("salinity","Seawater salinity [%]", None)
        self.AddParameter("charFile",".pkl file that hold the noise characterization",os.getenv('PONESRCDIR')+"/data/NoiseCharacterization.pkl")
        self.AddParameter("inputTreeName","Name of the input tree that holds physics data to span the noise generation over","I3Photon")
        self.AddParameter("outputTreeName","Name of the tree that the noise photons will be written to","K40")
        self.AddParameter("startPadding","Time padding before the first physics photon [ns]",1000)
        self.AddParameter("endPadding","Time padding after the last physics photon [ns]",10000)
        self.AddParameter("skipSingles","Choose whether to ignore events that only give 1 photon",False)
        self.AddParameter("pmtHitsOnly","Choose whether to only write photons that landed on a PMT",False)
        self.AddParameter("skipOneFold","Choose whether to ignore events that only hit 1 PMT. If this is set to True, pmtHitsOnly must also bet set to True.",False)



    def Configure(self):

        self.SEED        = self.GetParameter("SEED")
        self.salinity    = self.GetParameter("salinity")
        self.charFile    = self.GetParameter("charFile")
        self.inTree      = self.GetParameter("inputTreeName")
        self.outTree     = self.GetParameter("outputTreeName")
        self.startPad    = self.GetParameter("startPadding")
        self.endPad      = self.GetParameter("endPadding")
        self.skipSingles = self.GetParameter("skipSingles")
        self.pmtHitsOnly = self.GetParameter("pmtHitsOnly")
        self.skipOneFold = self.GetParameter("skipOneFold")

        #---------------------------------------------
        # Fixed parameters
        self.domR    = 0.2159

        if self.skipOneFold == True:
            if self.pmtHitsOnly == False:
                raise Exception("pmtHitsOnly must be set to True if skipOneFold is True")
            if self.skipSingles == False:
                self.skipSingles == True

        #---------------------------------------------
        # If only interested in PMT hits, set up the DOM geometry.
        # Here we build a cartesian matrix of the PMT positions
        # in the PONE DOM. The PMT indexing here should match
        # the PMT indexing used in PONEDOMLauncher.
        thetaVals = list(np.pi/180*np.array([65., 32.5, 115., 147.5]))
        phiVals1  = list(np.pi/180*np.array([360., 270., 180., 90.]))
        phiVals2  = list(np.pi/180*np.array([315., 225., 135., 45.]))

        # Theta and phi arrays of coordinates for each of the PMTs in the DOM
        thetaList = np.array([0. for i in np.arange(16)])
        phiList   = np.array([0. for i in np.arange(16)])

        pmtIndex = 0
        for t in np.arange(len(thetaVals)):
            even = True if t == 0 or t == 3 else False
            for p in np.arange(4):
                thetaList[pmtIndex] = thetaVals[t]
                if even:
                    phiList[pmtIndex] = phiVals1[p]
                else:
                    phiList[pmtIndex] = phiVals2[p]
                pmtIndex += 1

        # Convert theta and phi to Cartesian coordinates and normalize
        xList = np.multiply(np.sin(thetaList),np.cos(phiList))
        yList = np.multiply(np.sin(thetaList),np.sin(phiList))
        zList = np.cos(thetaList)

        self.pmtMatrix = np.array([xList,yList,zList]).T # Cartesian matrix
        del xList, yList, zList, thetaList, phiList, thetaVals, phiVals1, phiVals2



    # Function for checking whether a photon hit a PMT in the DOM
    #
    # @ Params:
    #           photonPosition - I3Position of the photon
    #
    # @ Return:
    #           Index of the PMT that was hit or None if no PMT was hit
    #
    def HitCheck(self, photonPosition):
        # Define additional PMT and module parameters
        pmtrad = (55.)*1E-3 # Effective PMT area is up to the size of the gel pad [m]
    
        angles = self.pmtMatrix.dot(np.array([photonPosition[0],photonPosition[1],photonPosition[2]])/np.linalg.norm(photonPosition))
        pmtHit = np.where(np.logical_and(self.domR*np.sin(np.arccos(angles))<=(pmtrad),angles>=0))[0]
    
        if len(pmtHit)==1: # One PMT hit, figure out which
            return pmtHit[0]
        elif len(pmtHit)>1: # Error, should not happen
            raise IndexError("Photon is hitting two (or more) PMTs at the same time!",pmtHit,photonPosition,angles[pmtHit])
        else: # No hit
            return None



    # Generates an array of K40 noise hits. This modules uses a characterization
    # stored in (NoiseGen.pkl). Be sure to check that the characterizaiton
    # is correct and corresponds to the correct DOM.
    # 
    # @Params:
    #           duration   - Time span over which to generate 40K events given in ns
    #           timeOffset - Offset from time = 0 by which to shift event times
    # 
    # @Return:
    #           Returns a numpy array of I3Photons
    def GenerateFrameNoise(self, duration, timeOffset):

        #---------------------------------------------
        # Load in the K40 Characterization
        k40 = unpickle(self.charFile)

        #---------------------------------------------
        # Generate noise and write to files

        # Set a seed
        if self.SEED:
            _=np.random.seed(self.SEED)

        # Calculate an event rate based on salinity if needed
        if self.salinity:
            eRate = getMeanEventRate(self.salinity)
        else:
            # If salinity is not specified just use the known rate
            # Hard coded here to save time
            eRate = 0.002540381
            #eRate = getMeanEventRate()

        # This is used as a Minor ID to track true coincidences
        decayNum = 1
        
        #---------------------------------------------
        # Get the timestamps of all the events in the frame
        frameTimes = k40.getFrameEvents(duration, eRate)
        frameTimes = frameTimes + timeOffset
        numEvents = len(frameTimes)
        
        #---------------------------------------------
        # Get the direction, distance to, and number of photons for each event
        eDirs = k40.getDir(num = numEvents)
        eDs = k40.getD(num =  numEvents)
        ePs = k40.getNumPs(eDs, num = numEvents)

        #---------------------------------------------
        # If skiping single photon events, remove these from the arrays created above
        if(self.skipSingles):
            singleIndices = np.where(ePs==1)[0]
            frameTimes = np.delete(frameTimes, singleIndices)
            eDirs = np.delete(eDirs, singleIndices)
            eDs = np.delete(eDs, singleIndices)
            ePs = np.delete(ePs, singleIndices)

            numEvents = len(frameTimes)

        # Initialize array that hold photons and masks to remove any photons that
        # were not part of a coincidence or did not hit a PMT
        photons = np.array([None for j in np.arange(np.sum(ePs))])
        photonNum = 0 # Just an index for adding photons to the np array

        hitPMTs = np.array([None for j in np.arange(np.sum(ePs))], dtype='object')
        
        #---------------------------------------------
        # Generate each detected decay
        for t in np.arange(numEvents):
            eventIndices = np.arange(photonNum, photonNum + ePs[t])
            
            #---------------------------------------------
            # Determine a location for the decay (Radius and Direction)
            dIndex, d = eDs[t]
            decayDirection = eDirs[t]
            
            #---------------------------------------------
            # Determine the number of photons that will reach the DOM
            numPs = ePs[t]
            
            #---------------------------------------------
            # Set up the rotation quaternion
            # We want to rotate from [0,0,1] to the decay direction
            zAxis = np.array([0,0,1])
            
            # Get the rotation axis vector by taking the cross product
            rAx = np.cross(decayDirection, zAxis)
            rAx = rAx/np.linalg.norm(rAx)
            
            # Get the angle we need to rotate by
            # Use the projections of the two vectors on the plane of rotation and
            # take the arccos of the dot product of the projections to get the angle
            
            startDirProj = zAxis - zAxis.dot(rAx) * rAx
            decayDirProj = decayDirection - decayDirection.dot(rAx) * rAx
            
            rotAng = 2*np.pi - np.arccos(startDirProj.dot(decayDirProj)) # Angle in radians

            # Define the quaternion object
            quat = myQuaternion(rAx, rotAng)
            
            if numPs < 2:
                # Choose an incident direction vector assuming the
                # hit is at [0,0,1]
                incDir = k40.getSingleInc()
                
                # Choose a wavelength
                wl = k40.getWl(d)
                
                # Rotate the hit position and direction by the quaternion
                hitPosition = quat.rotate(zAxis)
                hitDirection = quat.rotate(incDir)
                
                # Create the photon
                photons[photonNum] = (makePhoton(hitPosition, hitDirection, frameTimes[t], wl, decayNum))
                if self.pmtHitsOnly:
                    hitPMTs[photonNum] = self.HitCheck(hitPosition)
                photonNum += 1
                
                decayNum += 1
                continue
            
            #---------------------------------------------
            # Set up the projected line of hits
            
            # Get endpoints
            start, end = k40.getLinPoints(numPs)
            # Get line length
            linLen = np.linalg.norm(end-start)
            # Line longitudinal direction
            lDir = (end-start)/linLen
            # Line transverse direction
            tDir = np.array([-lDir[1], lDir[0]])
            
            offsetRIndex = k40.getOffsetRIndex(start, end)
            offsetDIndex = k40.getOffsetDIndex(d)
            
            # Get the transverse standard deviation
            tStdv = np.random.choice(k40.linStdvBinCentres, 1, p=k40.linSpreadDist[offsetDIndex])[0]
            # tStdv = k40.getStdv(offsetDIndex) # Removed to reduce function call overhead
            
            #---------------------------------------------
            # Loop over ever photon
            
            firstPhoton = True
            for p in np.arange(numPs):
                
                # Set up the hit position
                lPos = k40.smear(np.random.choice(k40.linLBinCentres, 1,
                                 p=k40.linPointDist[offsetRIndex])[0], k40.linLBinCentres)
                #lPos = k40.getLinPos(offsetRIndex)# Removed to reduce function call overhead
                tOffset = np.random.normal(0, tStdv)
        
                hitPos = start + lPos*linLen*lDir + tOffset*tDir
    
                # If the hit is out of the circle then just put it on the edge
                if np.linalg.norm(hitPos) > 1:
                    hitPos = 0.99*hitPos/np.linalg.norm(hitPos) # 0.99 so extending to 3D sqrt doesn't give err
        
                # Get the z coordinate based on x and y
                if hitPos[0]**2 - hitPos[1]**2 >=1.:
                    hitZ = 0.
                else:
                    hitZ = np.sqrt(1. - hitPos[0]**2 - hitPos[1]**2)
                
                tempHitPos = self.domR * np.array([hitPos[0], hitPos[1], hitZ])
                
                # Get the hit zenith based on the generated position
                zenith = (180/np.pi)*np.arccos(hitZ)
                
                # Get the number of scatters
                numScatters = np.random.choice(np.arange(4), 1, p=k40.sProbs[dIndex])[0]
                #k40.getS(dIndex) # Removed to reduce function call overhead
                
                # Choose a wavelength
                wl = k40.getWl(d)
                
                # Determine the incident direction
                devAng = k40.getDevAng(zenith, numScatters)
                incDir = k40.getIncDirNear(tempHitPos, d, devAng)
                
                # Rotate the hit position and direction by the quaternion
                hitPosition = quat.rotate(tempHitPos)
                hitDirection = quat.rotate(incDir)
                
                # Choose the hit time
                if firstPhoton:
                    hTime = frameTimes[t]
                    firstPhoton = False
                else:
                    hTime = k40.getNextTime(frameTimes[t], numScatters, timeOffset, timeOffset + duration)
                
                # Create a photon and add it to the list
                photons[photonNum] = (makePhoton(hitPosition, hitDirection, hTime, wl, decayNum))
                if self.pmtHitsOnly:
                    hitPMTs[photonNum] = self.HitCheck(hitPosition)
                photonNum += 1
            
            #---------------------------------------------
            # For each decay remove any events that give less than a two-fold coincidence
            if self.skipOneFold:
                eventPMTs = hitPMTs[eventIndices]
                eventPMTs = eventPMTs[eventPMTs != np.array(None)]
                if len(np.unique(eventPMTs)) < 2:
                    hitPMTs[eventIndices] = None

            decayNum += 1

        #---------------------------------------------
        # Now that we have all the photons as well as the mask of which we want to keep.
        # Apply the mask and return the remaining photons or just return all the photons
        # if we don't care about the PMTs
        if self.pmtHitsOnly:
            photonMask = np.where(hitPMTs != None, True, False)
            pmts = hitPMTs[np.where(hitPMTs != None)]
            return photons[photonMask], pmts
        else:
            return photons, None
    


    def DAQ(self, frame):

        photonDOMMap = frame[self.inTree]

        # These are used otherwise
        newPhotonMap = {}
        outputPulseMap = None
        outputname = ""
        if not self.pmtHitsOnly:
            outputPulseMap = simclasses.I3PhotonSeriesMap()
            outputname = '_Photons'
        else :
            outputPulseMap = dataclasses.I3RecoPulseSeriesMap()
            outputname = '_pmtsplit'

        for omkey in photonDOMMap.keys():
            # Determine the duration
            firstPhotonTime = photonDOMMap[omkey][0].time
            lastPhotonTime = photonDOMMap[omkey][-1].time

            # Add padding to the generation range. For insance the default is to add 1 us
            # to the beginning and 10 us to the end of the time span. -> 11000 in ns as padding
            timeSpan = (lastPhotonTime - firstPhotonTime) + self.startPad + self.endPad
            timeOffset = firstPhotonTime - self.startPad

            # Make a new object to store the noise hits
            
            newomkey = NoPMTKey(omkey)
            newPhotonMap[newomkey] = []

            noisePhotons, noisePMTs = self.GenerateFrameNoise(timeSpan, timeOffset)

            # If we don't care about PMTs just take the I3Photons and make an I3PhotonSeriesMap
            if not self.pmtHitsOnly:
                #noisePhotons = self.GenerateFrameNoise(10000, 0) # FOR EASY TESTING COMPARISON
                if len(noisePhotons) > 0:
                    outputPulseMap[omkey] = simclasses.I3PhotonSeries(noisePhotons)

            
            # If we do care about PMTs we need to turn the photons into I3RecoPulses and
            # index with the hit PMT
            else:
                for i in np.arange(len(noisePhotons)):
                    newPhotonMap[newomkey].append((noisePhotons[i].time, noisePMTs[i]))
                    pmtkey = AddPMTKey(newomkey, int(noisePMTs[i]))
                    if pmtkey not in outputPulseMap.keys():
                        outputPulseMap[pmtkey] = dataclasses.I3RecoPulseSeries()
                    rPulse = dataclasses.I3RecoPulse()
                    rPulse.time = noisePhotons[i].time
                    rPulse.charge = 1.0
                    outputPulseMap[pmtkey].append(rPulse)

                # Don't think I write this one to the I3 file. I get an error if I do. I think
                # we need the newPhotonMap to pass to the ApplyPMTResponse which will then give
                # TTS and pulse combining I think.
        frame[self.outTree+outputname] = outputPulseMap
        self.PushFrame(frame)

