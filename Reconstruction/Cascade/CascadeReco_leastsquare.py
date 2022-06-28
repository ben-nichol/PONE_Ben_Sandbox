#!/usr/bin/env python                                                                                                
# This is meant to be a slightly more robust approach to reconstruction of a muon event.                               
# The physics and likelihood model is heavily based off of the ICECUBE model and can be found at                      
# "https://publications.ub.uni-mainz.de/theses/volltexte/2014/3869/pdf/3869.pdf"                                     
# The time residuals are computed by myself though. The techniques used are detailed in a text document I have somewhere -dg
 
# Import some useful ICECUBE modules                                                                                  
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey 
from icecube.dataclasses import I3Particle 
import numpy as np
import time
from scipy import special as sp                        # For the Gamma function 
from scipy import optimize as op
from scipy import interpolate as interp
import sys
#from iminuit import Minuit
import argparse
import math as m
from Utilities.DOMUtility import AddPMTKey
from Utilities.OpticalParameters import c, ngroup

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,_minpos):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    distance = data
    geo_doms = domsUsed
    minpos = _minpos

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vx,vy,vz):

        vertex = dataclasses.I3Position(vx,vy,vz)
        leastsquare = 0.0
        deltal = np.sqrt((minpos.x-vertex.x)**2.0+(minpos.y-vertex.y)**2.0+(minpos.z-vertex.z)**2.0)
        for dom in distance.keys() :
            dompos = geo_doms[dom].position
            vertex_to_dom  = np.sqrt((dompos.x-vertex.x)**2.0+(dompos.y-vertex.y)**2.0+(dompos.z-vertex.z)**2.0)
            for dist in distance[dom] :
                leastsquare += dist[1]*(vertex_to_dom - dist[0] - deltal)**2.0

        return leastsquare

    return likelihoodFunction

def GetDistances(pulse_series,geo_doms):                                 

    c_n = c/ngroup                                     # light in water

    mintime = 10000.
    mintime_dom = None

    for domkey in pulse_series.keys() :
        for pulse in pulse_series[domkey] :
            if pulse.time < mintime :
                mintime = pulse.time
                mintime_dom = domkey

    T0 = mintime

    minpos = geo_doms[mintime_dom].position
    distance = {}
    for domkey in pulse_series.keys() :
        nopmtkey = AddPMTKey(domkey,1)
        if nopmtkey not in distance.keys() :
            distance[nopmtkey] = []
        for pulse in pulse_series[domkey] :
            dist = abs(pulse.time - mintime)*c_n
            distance[nopmtkey].append((dist,pulse.charge))

    return T0, minpos, distance

class CascadeReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo
        self.PushFrame(frame)

    def Physics(self,frame):

        if not frame.Has(self.pulseseries) :
            self.PushFrame(frame)
            return

        data = frame[self.pulseseries]

        T0, minpos, distance = GetDistances(data,self.domsUsed)
        vertex = minpos
        qFunctor = LikelihoodFunctor(distance,self.domsUsed,minpos)

        # Minimize using scipy
        def func(x):
            vx, vy, vz = x
            return qFunctor(vx,vy,vz)

        solution = op.minimize(fun=func, 
                               x0=np.array([vertex.x,vertex.y,vertex.z]), 
                               method='Nelder-Mead')

        q = dataclasses.I3Position(solution.x[0],solution.x[1],solution.x[2])

        totalcharge = 0.0
        totalcharge_dir = 0.0
        T0 = 0.0
        c_n = c/ngroup
        direc = [0.0,0.0,0.0]

        for dom in data.keys():
            dompos = self.domsUsed[dom].position
            Vertex_to_dom = np.sqrt((dompos.x-q.x)**2.0+(dompos.y-q.y)**2.0+(dompos.z-q.z)**2.0)
            if Vertex_to_dom > 0.0 :
                _direct = [(dompos.x-q.x)/Vertex_to_dom,(dompos.y-q.y)/Vertex_to_dom,(dompos.z-q.z)/Vertex_to_dom]
            for pulse in data[dom] :
                direc[0] += pulse.charge*_direct[0]
                direc[1] += pulse.charge*_direct[1]
                direc[2] += pulse.charge*_direct[2]
                T0 += pulse.charge*(pulse.time - Vertex_to_dom/c_n)
                totalcharge += pulse.charge
                if Vertex_to_dom > 0.0 :
                    totalcharge_dir += pulse.charge
        T0 /= totalcharge
        direc[0] /= totalcharge_dir
        direc[1] /= totalcharge_dir
        direc[2] /= totalcharge_dir
        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.Cascade

        # record on particle whether reconstruction was successful
        if solution.success == True:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
        recoParticle.speed = c
        recoParticle.pos = q
        recoParticle.dir = dataclasses.I3Direction(direc[0],direc[1],direc[2])
        recoParticle.time = T0 

        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_leastsquare"] =  dataclasses.I3Double(solution.fun)
        self.PushFrame(frame)    

