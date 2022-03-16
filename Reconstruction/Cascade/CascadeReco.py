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
from scipy import special as sp                        # For the Gamma function 
from scipy import optimize as op
from scipy import interpolate as interp
import sys
#from iminuit import Minuit
import argparse
import math as m
import random as rand
from Utilities.RecoUtility import GetPhotonTravelTime
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.OpticalParameters import c, n, ngroup, tau

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,_tau):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed

    vx = 0.0
    vy = 0.0
    vz = 0.0 
    v =np.array([0.0,0.0,0.0])

    thistau = _tau                                        # absorbtion length     

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vx,vy,vz,t0):
        dark = 1.e-16

        vertex = dataclasses.I3Position(vx,vy,vz)
        sum_nloglike = 0.0
        for dom in pulse_series.keys() :
            domkey =  dom #NoPMTKey(dom) fixed with GCD
            dc,t = GetPhotonTravelTime([geo_doms[domkey].position.x,geo_doms[domkey].position.y,geo_doms[domkey].position.z],[vertex.x,vertex.y,vertex.z])
            p_charge = np.exp(-dc/thistau)/max(dc,0.25)
            for pulse in pulse_series[dom] :
                charge = 1.0
                cpandel_out = pdf(pulse.time - t0 - t ,dc)
                if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
                    charge = pulse.charge
                sum_nloglike -= charge*np.log(cpandel_out*p_charge+dark)
                sum_nloglike -= charge*min(0.0,pulse.time - t0 - t)

        return sum_nloglike

    return likelihoodFunction

def GetVertexTime(pulse_series,geo_doms):                                 

	c_n = c/ngroup                                     # light in water
	ismc = False
	if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
		ismc = True
	
	totalcharge = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0

	for domkey in pulse_series.keys() :
		domkey_nopmt =  domkey #NoPMTKey(domkey) fixed with GCD
		for pulse in pulse_series[domkey] :
			totalcharge += pulse.charge
			vx += geo_doms[domkey_nopmt].position.x*pulse.charge
			vy += geo_doms[domkey_nopmt].position.y*pulse.charge
			vz += geo_doms[domkey_nopmt].position.z*pulse.charge

	if totalcharge < 5.0 :
		return 0.0, dataclasses.I3Position(0.0,0.0,0.0), totalcharge
	vertex = dataclasses.I3Position(vx/totalcharge,vy/totalcharge,vz/totalcharge)

	T0 = 0.0

	for domkey in pulse_series.keys() :
		domkey_nopmt =  domkey#NoPMTKey(domkey)fixed with GCD
		for pulse in pulse_series[domkey] :
			dx = vertex.x - geo_doms[domkey_nopmt].position.x
			dy = vertex.y - geo_doms[domkey_nopmt].position.y
			dz = vertex.z - geo_doms[domkey_nopmt].position.z
			dist = np.sqrt(dx*dx+dy*dy+dz*dz)
			T0 += pulse.time - dist/c_n
            
	if totalcharge < 5.0 :
		return T0, vertex, totalcharge
	T0 /= totalcharge
	T0 -= 5.0
	return T0, vertex, totalcharge

class CascadeReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("tau","optical attenuation length.",tau)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")

        self.tau = self.GetParameter("tau")

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo
        self.PushFrame(frame)

    def DAQ(self,frame):

        if not frame.Has(self.pulseseries) :
            self.PushFrame(frame)
            return

        data = frame[self.pulseseries]

        qFunctor = LikelihoodFunctor(data,self.domsUsed,self.tau)
        T0, vertex, totalcharge = GetVertexTime(data,self.domsUsed)

        if totalcharge < 5.0 :
            return 

        # Minimize using scipy
        def func(x):
            vx, vy, vz,t0 = x
            return qFunctor(vx,vy,vz,t0)

        solution = op.minimize(fun=func, 
                               x0=np.array([vertex.x,vertex.y,vertex.z,T0]), 
                               method='Nelder-Mead')

        q = dataclasses.I3Position(solution.x[0],solution.x[1],solution.x[2])

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
        recoParticle.time = solution.x[3]

        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(solution.fun)
        self.PushFrame(frame)    

