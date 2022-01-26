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
from Utilities.PandelPDFs import cpandel as pdf               # This module is used to store the pdf
from scipy import optimize as op
import sys
import argparse
import math as m
from Utilities.RecoUtility import GetGeoTime
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.OpticalParameters import c, n, ngroup

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,vertexrad):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed

    vx = 0.0
    vy = 0.0
    vz = 0.0 
    v =np.array([0.0,0.0,0.0])

    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 18.949132224466762                                        # time parameter that has to be fit using simulations or data      
    vertexRad = vertexrad
    # min time index for the first hit PMT

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vtheta, vphi, theta, phi, t0):
        dark = 1.e-8

        vertex = dataclasses.I3Position(vertexRad*np.sin(vtheta)*np.cos(vphi),vertexRad*np.sin(vtheta)*np.sin(vphi),vertexRad*np.cos(vtheta))
        direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))

        sum_nloglike = 0.0
        for dom in pulse_series.keys() :
            domkey =  NoPMTKey(dom) 
            d,dc,t = GetGeoTime([geo_doms[domkey].position.x,geo_doms[domkey].position.y,geo_doms[domkey].position.z],
                                [vertex.x,vertex.y,vertex.z],
                                [direction.x,direction.y,direction.z])
            p_charge = np.exp(-d/tau)/max(dc,0.25)
            for pulse in pulse_series[dom] :
                charge = 1.0
                time_r = pulse.time - t0 - t
                cpandel_out = pdf(time_r ,d)
                if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
                    charge = pulse.charge                
                if time_r < 0 :
                    sum_nloglike -= charge*np.log(cpandel_out*p_charge+dark) + time_r
                else :
                    sum_nloglike -= charge*np.log(cpandel_out*p_charge+dark)

        return sum_nloglike

    return likelihoodFunction

def GetVertexTime(vertex,direction,pulse_series,geo_doms):                                 

    totalcharge = 0.0
    MaxChargeDOM = None
    maxCharge=0.0
    DOMCharge = {}

    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n) 
    for dom in pulse_series.keys() :
        totalcharge = 0.0
        for pulse in pulse_series[dom] :
            totalcharge += pulse.charge
        if NoPMTKey(dom) in DOMCharge.keys() :
            DOMCharge[NoPMTKey(dom)] += totalcharge
        else :
            DOMCharge[NoPMTKey(dom)] = totalcharge

    #time of largest pulse
    maxCharge=0.0
    maxCharge_time = 0.0
    for dom in DOMCharge.keys() :
        if DOMCharge[dom] > maxCharge :
            maxCharge = DOMCharge[dom]
            MaxChargeDOM = dom

    if type(MaxChargeDOM) != type(OMKey(0,0,0)) :
        return 7200

    DOMPos = geo_doms[MaxChargeDOM].position

    maxCharge=0.0
    maxCharge_time = 0.0
    for domkey in pulse_series.keys():
        if (domkey.string == MaxChargeDOM.string) and (domkey.om == MaxChargeDOM.om):
            for pulse in pulse_series[domkey] :
                if pulse.charge > maxCharge :
                    maxCharge = pulse.charge
                    maxCharge_time = pulse.time

    x = DOMPos.x - vertex.x
    y = DOMPos.y - vertex.y
    z = DOMPos.z - vertex.z
    # Compute (\vec{r} - vec{x}) dot \vec{v}
    dotprod = x*direction.x + y*direction.y + z*direction.z
    # Compute the final vector components
    # Compute t_i,c and d_i,c
    dc = np.sqrt(x*x + y*y + z*z - dotprod*dotprod)
    #time to travel to closest approach
    tc = dotprod/c

    # Now we find the time of the photon emission
    _tc = tc - dc/(np.tan(theta_c)*c)
    # The first component of the geometric time
    d = dc/np.sin(theta_c)
    t_geo = d/c_n   
    # The total geometric time
    t_geo = t_geo + _tc
    # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
    return maxCharge_time - t_geo

class TrackReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("seedtrack","Track to seed fit","linefit")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("vertexRad","Radius to put vertex at",550.)
        self.AddParameter("UseMC","Use MC Truth Track to seed",False)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.useMC = self.GetParameter("UseMC")
        self.domsUsed = {}
        # Some quantities that are environment dependent
        self.theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 18.949132224466762                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo

        maxradius = 0.0
        for dom in self.domsUsed :
            pos = self.domsUsed[dom].position
            radius = np.sqrt(pos.x**2.0+pos.y**2.0+pos.z**2.0)
            maxradius = max(maxradius,radius)

        self.vertexRad = maxradius + 100.0

    def DAQ(self,frame): 
        data = frame[self.pulseseries]

        linefit = frame[self.seedtrack]

        direction = dataclasses.I3Direction(linefit.dir.x,linefit.dir.y,linefit.dir.z) 

        qFunctor = LikelihoodFunctor(data,self.domsUsed,self.vertexRad)

        p_2 = linefit.pos.x**2.0+linefit.pos.y**2.0+linefit.pos.z**2.0
        pd = (linefit.pos.x*direction.x+linefit.pos.y*direction.y+linefit.pos.z*direction.z)
        r_2 = self.vertexRad**2.0

        if pd**2.0-p_2+r_2 < 0.0 :
            return

        L = -pd - np.sqrt(pd**2.0-p_2+r_2)

        vertex = dataclasses.I3Position(linefit.pos.x+L*direction.x,linefit.pos.y+L*direction.y,linefit.pos.z+L*direction.z)
        
        VTheta = vertex.theta
        VPhi = vertex.phi

        T0 = GetVertexTime(vertex,direction,data,self.domsUsed)

        # Minimize using scipy
        def func(x):
            vtheta, vphi, theta, phi, t0 = x
            return qFunctor(vtheta, vphi, theta, phi, t0)
        solution = op.minimize(fun=func, 
                               x0=np.array([VTheta, VPhi, direction.theta, direction.phi, T0]), 
                               method='Nelder-Mead')

        vx = self.vertexRad*np.sin(solution.x[0])*np.cos(solution.x[1])
        vy = self.vertexRad*np.sin(solution.x[0])*np.sin(solution.x[1])
        vz = self.vertexRad*np.cos(solution.x[0])
        q = dataclasses.I3Position(vx,vy,vz)
        phi = solution.x[3] 
        theta = solution.x[2]
        u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.InfiniteTrack
                
        # record on particle whether reconstruction was successful
        #if minimizer.get_fmin()["is_valid"]:
        if solution.success == True:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
        recoParticle.dir = u
        recoParticle.speed = c
        recoParticle.pos = q
        recoParticle.time = solution.x[4]

        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(solution.fun)
        frame[self.output+"_seed_llhval"] = dataclasses.I3Double(qFunctor(VTheta, VPhi, direction.theta, direction.phi, T0))
        
        self.PushFrame(frame)    

