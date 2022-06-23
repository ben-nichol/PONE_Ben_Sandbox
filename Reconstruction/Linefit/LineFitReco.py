#!/usr/bin/env python

'''
SimAnalysis contains helper functions commonly used when writing
scripts to analyze results from simulation sets. Mainly, it allows
to easily decompose geometries into smaller sets by specifying a list
of omkeys, and allows to impose harder cuts on simulations so that different
geometries/cuts could be tried out without running the simulation again.
'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import numpy as np
from numpy import linalg as la
from icecube.phys_services import I3Calculator
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
class LineFitReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("inputseries","Input pulse series","MCPESeriesMap")
        self.AddParameter("output","Output track name.","linefit")
        self.AddParameter("hitThresh","Threshold for number of pulses in DOM",1)
        self.AddParameter("domThresh","Threshold for number of DOMs",6)
        self.AddParameter("vertRadius","radius to readust the vertext to",200.)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.hitThresh = self.GetParameter("hitThresh")
        self.domThresh = self.GetParameter("domThresh")
        self.input = self.GetParameter("inputseries")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertRadius")

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo

        maxradius = 0.0
        for dom in self.domsUsed.keys() :
            pos = self.domsUsed[dom].position
            radius = np.sqrt(pos.x**2.0+pos.y**2.0+pos.z**2.0)
            maxradius = max(maxradius,radius)

        self.vertexRad += maxradius 

        self.PushFrame(frame)

    # A function that determines whether a frame passes the cut or not. The
    # function checks whether a threshold number of DOMs passed given a 
    # threshold number of hits needed in the 20ns window and a threshold
    # residual to count the hits 
    # 
    # @Param:
    # frame:            The frame in question   
    # domsUsed:         List of omkeys used for the analysis. Allows to look 
    #                   at smaller geometries from a larger geometry sim file
    # hitThresh:        Hit threshold for a single DOM
    # domThresh:        Number of passed DOMs needed to pass the frame 
    # maxResidual:      Maximum time residual allowed for a DOM to be considered
    # GeoMap:           An I3OMGeoMap object that maps the omkeys in the
    #                   geometries to their I3OMGeo objects 
    # 
    # @Return: 
    # A boolean variable indicating whether the frame passed or not
    def passFrame(self,frame):
        mcpeMap = frame[self.input]
    
        domCount = 0
        for dom in mcpeMap.keys():
            if len(mcpeMap[dom]) >= self.hitThresh:
                domCount += 1
        
        if domCount >= self.domThresh:
            return True
        #print("too few doms")
        return False


    # parses the hit information requires for the linefit reconstruction and returns
    # a list of tuples containing the required hit information
    # NOTE: function assumes the frame already contains the significant MCPESeriesMap
    #       and tries to call on the "MCPESeriesMap_significant_hits" key. If the key
    #       does not exist, the function will raise a ValueError. It is much more  
    #       efficient to have the significant hits MCPESeriesMap to use rather than
    #       having to search for the significant hits every time, so most methods rely
    #       on the frame key being there.   
    #
    # @Param:
    # frame:            The frame containing the event information   
    # geometry:         An I3Geometry object from the gcd file used to produce the
    #                   simulation data 
    # hitThresh:        Hit threshold for a single DOM
    #
    # @Return: 
    # A list of tuples. Each tuple contains information for a single hit (x,y,z,t)
    def getLinefitDataPoints(self,frame):
        if not frame.Has(self.input):
            raise ValueError("Frame does not contain " + self.input)
        mcpeMap = frame[self.input]
        data = []
        for omkey, mcpeList in mcpeMap:
            timeList = [mcpe.time for mcpe in mcpeList]
            npeList = [mcpe.charge for mcpe in mcpeList]
            if len(timeList)<1 or len(npeList)<1 :
              continue
            time = min(timeList)
            charge = sum(npeList)
            key = omkey #NoPMTKey(omkey)fixed with GCD
            position = self.domsUsed[omkey].position
            for i in range(len(timeList)):
                data.append( (position.x, position.y, position.z, time, charge) )
    
        return data

    # Given a list of x and y points, computes the parameters of a least squares fit
    # line
    #  
    # @Param:
    # x:                The x component of points in the fit   
    # y:                The y component of points in the fit
    # @Return: 
    # A tuple containing the slope and y intercept of least squares fit line
    def fitLeastSquaresLine(self,x, y):
        # for a given system X*c = y 
        # where X is a truncated Vandermond matrix (truncated in this case to be degree 1)
        #       c is the vector containing the polyniomial coefficients
        #       y is the vector containing all y values
        # The least squares solution is given as c = (X.T * X)^(-1) * (X.T * y)
        xMatrix = np.column_stack( (x, np.ones(len(x))) )
        yVector = np.array(y).T
        leastSquaresMatrix = np.matmul(xMatrix.T, xMatrix)
        #if la.det(leastSquaresMatrix) == 0:
        #    print leastSquaresMatrix, xMatrix
        leastSquaresVector = np.matmul(xMatrix.T, yVector)
        fitCoefficients = np.matmul( la.inv(leastSquaresMatrix), leastSquaresVector)

        slope = fitCoefficients[0]
        intercept = fitCoefficients[1]

        return slope, intercept

    # Given datapoints containing the hit informations (from getLineFitDatapoints),
    # the function computes the parameters for the linefit reconstruction of the 
    # particle 
    # 
    # @Param:
    # datapoints:       A list of tuples. Each tuple contains hit information in the
    #                   format (x,y,z,t,charge) 
    # 
    # @Return: 
    # A tuple containing the reconstructed particle's information. This is in the
    # of an I3Direction object for the particle's direction, a double for the particle's
    # speed, and an I3Position object for the particle vertex (position at t=0)  
    def Physics(self,frame):
        if not self.passFrame(frame) :
            #self.PushFrame(frame)
            return

        datapoints = self.getLinefitDataPoints(frame)
        x = [data[0] for data in datapoints]
        y = [data[1] for data in datapoints]
        z = [data[2] for data in datapoints]
        t = [data[3] for data in datapoints]
        charge = [data[4] for data in datapoints]

        weighted_time = np.sum(np.array(t) * np.array(charge))/np.sum(np.array(charge))

        xVelocity, x = self.fitLeastSquaresLine(t, x)
        yVelocity, y = self.fitLeastSquaresLine(t, y)
        zVelocity, z = self.fitLeastSquaresLine(t, z)

        direction = dataclasses.I3Direction(xVelocity, yVelocity, zVelocity)
        speed = np.sqrt(xVelocity**2 + yVelocity**2 + zVelocity**2)
        vertex = dataclasses.I3Position(x,y,z)

        #reset vertex to be at set radius, same as track fitter.

        radius_0 = vertex.x**2.0+vertex.y**2.0+vertex.z**2.0
        b = (vertex.x*direction.x+vertex.y*direction.y+vertex.z*direction.z)
        c = radius_0-self.vertexRad**2.0
        b2_4ac = b*b-c
        l = 0.0
        if b2_4ac > 0.0 :
            l = -b -np.sqrt(b2_4ac)
        else :
            l = -b
        
        vertex = dataclasses.I3Position(x+l*direction.x,y+l*direction.y,z+l*direction.z)

        linefit = dataclasses.I3Particle()
        linefit.shape = dataclasses.I3Particle.InfiniteTrack
        linefit.dir = direction                                                                                                                  
        linefit.speed = speed
        linefit.pos = vertex        
        linefit.time = weighted_time

        frame[self.output] = linefit
        self.PushFrame(frame)  


