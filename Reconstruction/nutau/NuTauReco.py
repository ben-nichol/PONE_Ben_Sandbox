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
from Reconstruction.llh.reco_pdfs import cpandel as pdf               # This module is used to store the pdf
from scipy import special as sp                        # For the Gamma function 
from scipy import optimize as op
import sys
#from iminuit import Minuit
import argparse
import math as m
import random as rand

# Geometric Time computation:
def GetGeoTime(position,vert) :
    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    x = position.x - vert.x
    y = position.y - vert.y
    z = position.z - vert.z
    dc = np.sqrt(x*x + y*y + z*z)
    t = dc/c_n
    return dc,t

def anisotropy(position,vert,direction):

    x1 = position.x - vert.x
    y1 = position.y - vert.y
    z1 = position.z - vert.z
    r1= m.sqrt(x1*x1+y1*y1+z1*z1)

    x2= direction.x
    y2= direction.y
    z2= direction.z
    zeta = 0.0
    if r1 > 0.0 :
    	zeta=np.arccos((x1*x2+y1*y2+z1*z2)/r1)
    

    return weight

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed

    vx = 0.0
    vy = 0.0
    vz = 0.0 
    v =np.array([0.0,0.0,0.0])

    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 18.949132224466762                                        # time parameter that has to be fit using simulations or data      

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vx,vy,vz,theta,phi,t0):
        dark = 1.e-16

        vertex = dataclasses.I3Position(vx,vy,vz)
        direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)) 
        sum_nloglike = 0.0
        for dom in pulse_series.keys() :
            domkey =  OMKey(dom.string, dom.om, 0) 
            dc,t = GetGeoTime(geo_doms[domkey].position,vertex)
            p_charge = np.exp(-dc/tau)/max(dc,0.25)
            anisotropyweight = anisotropy(geo_doms[domkey].position,vertex,direction)
            for pulse in pulse_series[dom] :
                charge = 1.0
                cpandel_out = pdf(pulse.time - t0 - t ,dc)
                if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
                    charge = pulse.charge
                sum_nloglike -= charge*np.log(cpandel_out*p_charge*anisotropyweight+dark)
                sum_nloglike -= charge*min(0.0,pulse.time - t0 - t)

        return sum_nloglike

    return likelihoodFunction

def GetVertexTime(pulse_series,geo_doms):                                 

	c = 0.299792458                                 # speed of light 
	n = 1.34
	ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
	c_n = c/ngroup                                     # light in water
	ismc = False
	if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
		ismc = True
	
	totalcharge = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0

	for domkey in pulse_series.keys() :
		domkey_nopmt =  OMKey(domkey.string, domkey.om, 0)
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
		domkey_nopmt =  OMKey(domkey.string, domkey.om, 0)
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

class NuTauReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","nutaufit")
        self.AddParameter("electronfile","","")
        self.AddParameter("taufile","","")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34  
        self.ngroup = 1.3555714017                                      # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.ngroup                                       # light in water
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 18.949132224466762                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data


        self.electrontable_x = list()
        self.electrontable_y = list()
        self.tautable_x = list()
        self.tautable_y = list()

        infile = open(self.GetParameter("electronfile"),"r")
        lines = infile.readlines()
        linecount = 0
        for line in lines :
            splitline = line.split(" ",100)
            self.electrontable_x.append(splitline[0])
            self.electrontable_y.append(splitline[1])
        infile.close()

        infile = open(self.GetParameter("taufile"),"r")
        lines = infile.readlines()
        linecount = 0
        for line in lines :
            splitline = line.split(" ",100)
            self.tautable_x.append(splitline[0])
            self.tautable_y.append(splitline[1])
        infile.close()

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def DAQ(self,frame): 
        data = frame[self.pulseseries]

        domsUsed = frame['I3Geometry'].omgeo

        qFunctor = LikelihoodFunctor(data,domsUsed)

        T0, vertex, totalcharge = GetVertexTime(data,domsUsed)

        if totalcharge < 5.0 :
            return 

        # Minimize using scipy
        def func(x):
            vx, vy, vz,theta,phi,t0 = x
            return qFunctor(vx,vy,vz,theta,phi,t0)
        solution = op.minimize(fun=func, 
                               x0=np.array([vertex.x,vertex.y,vertex.z,0.0,0.0,T0]), 
                               method='Nelder-Mead')

        vx = solution.x[0]
        vy = solution.x[1]
        vz = solution.x[2]
        q = dataclasses.I3Position(vx,vy,vz)
        theta = solution.x[3]
        phi = solution.x[4] 
        d = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))  
        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.Cascade
                
        # record on particle whether reconstruction was successful
        if solution.success == True:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
        recoParticle.dir = d
        recoParticle.speed = self.c
        recoParticle.pos = q
        recoParticle.time = solution.x[5]

        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(solution.fun)
        
        self.PushFrame(frame)    

