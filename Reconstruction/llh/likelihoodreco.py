#!/usr/bin/env python                                                                                                
# This is meant to be a slightly more robust approach to reconstruction of a muon event.                               
# The physics and likelihood model is heavily based off of the ICECUBE model and can be found at                      
# "https://publications.ub.uni-mainz.de/theses/volltexte/2014/3869/pdf/3869.pdf"                                     
# The time residuals are computed by myself though. The techniques used are detailed in a text document I have somewhere -dg
 
# Import some useful ICECUBE modules                                                                                  
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
from Reconstruction.llh.reco_pdfs import log_cpandel as pdf               # This module is used to store the pdf
from scipy import special as sp                        # For the Gamma function 
import sys
from iminuit import Minuit
import argparse
import math as m
from Reconstruction.llh.ChargeLikelihood import nLogLikelihood

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(self,data,domsUsed,vertexrad,prnt = False):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    pmt = []
    time = []
    charge = []

    for dom in pulse_series.keys() :
        for pulse in pulse_series[dom] :
            pmt.append([])
            pmt[-1].append(geo_doms[dom].position.x)
            pmt[-1].append(geo_doms[dom].position.y)
            pmt[-1].append(geo_doms[dom].position.z)
            time.append(pulse.time)
            charge.append(pulse.charge)

    c = 0.299792458                                 # speed of light 
    n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/n                                       # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 557                                       # time parameter that has to be fit using simulations or data      
    vertexRad = vertexrad
    # min time index for the first hit PMT
    min_index = np.argmin(time)
    
    # The computations from here on require we find the time and distance of closest approach, d_i,c and t_i,c
    def closestApproach(vtheta, vphi, theta, phi):
        # Compute vec{r} - vec{x}
        vx = vertexRad*np.cos(vphi)*np.sin(vtheta)
        vy = vertexRad*np.sin(vphi)*np.sin(vtheta)
        vz = vertexRad*np.cos(vtheta)
        x = pmt[:,0] - vx
        y = pmt[:,1] - vy
        z = pmt[:,2] - vz
        # Compute (\vec{r} - vec{x}) dot \vec{v}
        v = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        dotprod = x*v[0] + y*v[1] + z*v[2]
        # Compute the final vector components
        x = x - dotprod*v[0]
        y = y - dotprod*v[1]
        z = z - dotprod*v[2]
        # Compute t_i,c and d_i,c
        dc = np.sqrt(x*x + y*y + z*z)
        tc = dotprod/c
        return dc, tc

    def GetVertexTime(vtheta,vphi) :                                      
        x = pmt[:,0]                                                                  
        y = pmt[:,1]                                                                  
        z = pmt[:,2]                                                                  
                                                                                 
        mean_t = 0.0                                                                  
        mean_x = 0.0                                                                  
        mean_y = 0.0                                                                  
        mean_z = 0.0                                                                  
        sum_charge = 0.0;                                                             
                                                                                  
        for i in range(len(time)):                                                    
            mean_t = mean_t + charge[i]*time[i]                                         
            mean_x = mean_x + charge[i]*x[i]                                            
            mean_y = mean_y + charge[i]*y[i]                                            
            mean_z = mean_z + charge[i]*z[i]                                            
            sum_charge = sum_charge + charge[i]                                         
                                                                                         
            mean_t = mean_t/sum_charge                                                    
            mean_x = mean_x/sum_charge                                                    
            mean_y = mean_y/sum_charge                                                    
            mean_z = mean_z/sum_charge                                                    
                            
            vx = vertexRad*np.cos(vphi)*np.sin(vtheta)                                
            vy = vertexRad*np.sin(vphi)*np.sin(vtheta)                                
            vz = vertexRad*np.cos(vtheta)  

        dist_phot = 0.0                                                               
                                                                                         
        for i in range(len(time)):                                                    
            dist_phot = dist_phot + charge[i]*np.sqrt((mean_x-x[i])**2.0+(mean_y-y[i])**2.0+(mean_z-z[i])**2.0)                                                          
        dist_phot = dist_phot/sum_charge                                              
                                                                                          
        dist = np.sqrt((vx-mean_x)**2.0+(vy-mean_y)**2.0+(vz-mean_z)**2.0)            
        vertextime = mean_t - dist/0.3 - dist_phot/(0.3/1.35)                         
        return vertextime 
        
    # Given the linefit and other parameters 
    def computeResiduals(dc, tc, t0):
        # Now we find the time of the photon emission
        tc = tc - dc/(np.tan(theta_c)*c)
        # The first component of the geometric time
        d = dc/np.sin(theta_c) 
        t_geo = d/c_n
        # Apply our offset time to find the "true" closest approach time array. Multiply by 1E9 to change to nanoseconds
        tc = tc + t0
        # The total geometric time
        t_geo = t_geo + tc 
        # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
        t = time - t_geo
        return d, t

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vtheta, vphi, theta, phi, t0): 
        dc, tc = closestApproach(vtheta, vphi, theta, phi)
        d, t = computeResiduals(dc, tc, t0)
        vx = vertexRad*np.sin(vtheta)*np.cos(vphi)
        vy = vertexRad*np.sin(vtheta)*np.sin(vphi)
        vz = vertexRad*np.cos(vtheta)
        charge_out = nLogLikelihood(pmt,charge,vx,vy,vz,theta,phi)
        out = pdf(t,d)
        return np.sum(out) + charge_out

    return likelihoodFunction

class likelihoodreco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("GCDFile","GCD file.")
        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("seedtrack","Track to seed fit","linefit")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("vertexRad","Radius to put vertex at",500.)

        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.gcdfile = self.GetParameter("GCDFile")
        self.geometry = self.gcdfile.pop_frame()["I3Geometry"]
        self.domsUsed = self.geometry.omgeo
        self.vertexRad = self.GetParameter("vertexRad")

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = c/n                                       # light in water
        self.theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 15.                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data      

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 

        data = frame[self.pulseseries]
        # Clean the data to get rid of repeated events
        #data = clean_data(data)
        linefit = frame[self.seedtrack]
          
        qFunctor = LikelihoodFunctor(data,self.domsUsed,self.vertexRad) 
        vr = np.sqrt(linefit.pos.x**2.+linefit.pos.y**2.0+linefit.pos.z**2.0)
        VTheta = np.arccos(linefit.pos.z/vr)
        VPhi = 0.0
        if np.sin(VTheta) != 0.0 :
            VPhi = np.arccos(linefit.pos.x/(vr*np.sin(VTheta)))   
        T0 = qFunctor.GetVertexTime(VTheta,VPhi)
          
        minimizer = Minuit(qFunctor, 
                        t0=T0,
                        error_t0=1.0,
                        vtheta=VTheta,
                        error_vtheta=1.0,
                        limit_vtheta=(0.0,np.pi),
                        vphi=VPhi,
                        error_vphi=1.0,
                        limit_vphi=(0.0,2.0*np.pi),
                        phi=linefit.dir.phi,  
                        error_phi=1.0,
                        limit_phi=(0.0,2.0*np.pi),
                        theta=linefit.dir.theta,
                        error_theta=1.0,
                        limit_theta=(0.0,np.pi),
                        errordef=0.5,
                        )

        minimizer.migrad()

        solution = minimizer.values
            
        # For likelihood
        vx = self.vertexRad*np.sin(solution['vtheta'])*np.cos(solution['vphi'])
        vy = self.vertexRad*np.sin(solution['vtheta'])*np.sin(solution['vphi'])
        vz = self.vertexRad*np.cos(solution['vtheta'])
        q = dataclasses.I3Position(vx,vy,vz)
        phi = solution['phi'] 
        theta = solution['theta']
        u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.InfiniteTrack
                
        # record on particle whether reconstruction was successful
        if minimizer.get_fmin()["is_valid"]:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
        recoParticle.dir = u
        recoParticle.speed = c
        recoParticle.pos = q
        recoParticle.time = 0
            
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  minimizer.fval
        self.push_frame(frame)    

