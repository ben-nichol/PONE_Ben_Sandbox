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

# A function that returns a unit vector direction such that it forms the angle alpha with the given direction. 
# Works for alpha in (0, pi/2) exclusive
def rand_dir(theta_dir, phi_dir, alpha):
    theta = 2*np.pi*rand.uniform(0,1)
    ct_dir = np.cos(theta_dir)
    st_dir = np.sin(theta_dir)
    cp_dir = np.cos(phi_dir)
    sp_dir = np.sin(phi_dir)
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    ct = np.cos(theta)
    st = np.sin(theta)
    dir_x = cp_dir*ct_dir*sa*ct + st_dir*ca*cp_dir - sp_dir*sa*st
    dir_y = sp_dir*ct_dir*sa*ct + st_dir*ca*sp_dir + cp_dir*sa*st
    dir_z = ca*ct_dir - sa*ct*st_dir
    return np.array([dir_x, dir_y, dir_z])

# Geometric Time computation:
def GetGeoTime(position,vert,direction) :
    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n)
    x = position.x - vert.x
    y = position.y - vert.y
    z = position.z - vert.z
    dotprod = x*direction.x + y*direction.y + z*direction.z
    dc = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
    d = dc/np.sin(theta_c)
    t = d/c_n + dotprod/c - dc/(np.tan(theta_c)*c)
    return d,dc,t


# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,vertexrad):
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
            domkey =  OMKey(dom.string, dom.om, 0) 
            d,dc,t = GetGeoTime(geo_doms[domkey].position,vertex,direction)
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

def GetStartStop(vertex,direction,pulse_series,geo_doms) :
        n = 1.34
        theta_c = np.arccos(1./n)  

        start = 99999999999
        stop = 0.0

        for dom in pulse_series.keys() :
            key = OMKey(dom.string,dom.om,0)
            dompos = geo_doms[key].position

            x = dompos.x - vertex.x
            y = dompos.y - vertex.y
            z = dompos.z - vertex.z
            # Compute (\vec{r} - vec{x}) dot \vec{v}
            dotprod = x*direction.x + y*direction.y + z*direction.z
            # Compute the final vector components
            # Compute t_i,c and d_i,c
            dc = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
            l = dotprod - dc/np.tan(theta_c)

            if l < start :
                start = l
            if l > stop :
                stop = l

        return start, stop

def GetVertexTime(vertex,direction,pulse_series,geo_doms):                                 

    totalcharge = 0.0
    MaxChargeDOM = None
    maxCharge=0.0

    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n) 
    if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
        for dom in pulse_series.keys() :
            totalcharge = 0.0
            for pulse in pulse_series[dom] :
                totalcharge += pulse.charge
            if totalcharge > maxCharge :
                maxCharge = totalcharge
                MaxChargeDOM = OMKey(dom.string,dom.om,0)
    else :
        for dom in pulse_series.keys() :
            totalcharge = 0.0
            for pulse in pulse_series[dom] :
                totalcharge += 1.0
            if totalcharge > maxCharge :
                maxCharge = totalcharge
                MaxChargeDOM = OMKey(dom.string,dom.om,0)

    #time of largest pulse
    maxCharge=0.0
    maxCharge_time = 0.0
    if type(MaxChargeDOM) != type(OMKey(0,0,0)) :
        return 7200
    DOMPos = geo_doms[MaxChargeDOM].position
	
    for pulse in pulse_series[MaxChargeDOM] :
        #print("charge = "+str(pulse.charge)+" time = "+str(pulse.time))
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
    #print("tc = "+str(tc))

    # Now we find the time of the photon emission
    _tc = tc - dc/(np.tan(theta_c)*c)
    # The first component of the geometric time
    d = dc/np.sin(theta_c)
    t_geo = d/c_n   
    #print("t_geo = " + str(t_geo))
    # The total geometric time
    t_geo = t_geo + _tc
    #print("maxcharge_time = "+str(maxCharge_time))
    # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
    return maxCharge_time - t_geo

class likelihoodreco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("seedtrack","Track to seed fit","linefit")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("vertexRad","Radius to put vertex at",550.)
        self.AddParameter("UseMC","Use MC Truth Track to seed",False)
        self.AddParameter("group","Used to fix at true direction, vertex or time",False)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.useMC = self.GetParameter("UseMC")
        self.group = self.GetParameter("group")

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34  
        self.ngroup = 1.3555714017                                      # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.ngroup                                       # light in water
        self.theta_c = np.arccos(1./self.n)                       # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 18.949132224466762                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def DAQ(self,frame): 
        tic = time.perf_counter() #Check start time
        data = frame[self.pulseseries]
        # Clean the data to get rid of repeated events
        #data = clean_data(data)

        linefit = frame[self.seedtrack]

        domsUsed = frame['I3Geometry'].omgeo

        direction = dataclasses.I3Direction(linefit.dir.x,linefit.dir.y,linefit.dir.z)
        
        #linefit = frame[self.seedtrack]
        #domsUsed = frame['I3Geometry'].omgeo 

        qFunctor = LikelihoodFunctor(data,domsUsed,self.vertexRad)

        p_2 = linefit.pos.x**2.0+linefit.pos.y**2.0+linefit.pos.z**2.0
        pd = (linefit.pos.x*direction.x+linefit.pos.y*direction.y+linefit.pos.z*direction.z)
        r_2 = self.vertexRad**2.0

        if pd**2.0-p_2+r_2 < 0.0 :
            return

        L = -pd - np.sqrt(pd**2.0-p_2+r_2)

        vertex = dataclasses.I3Position(linefit.pos.x+L*direction.x,linefit.pos.y+L*direction.y,linefit.pos.z+L*direction.z)
        # Load muon information
        MMCTrackList = frame['MMCTrackList']
        muon = MMCTrackList[0].GetI3Particle()

        if self.useMC :
            MMCTrackList = frame['MMCTrackList']
            linefit = MMCTrackList[0].GetI3Particle()
            direction = dataclasses.I3Direction(linefit.dir.x,linefit.dir.y,linefit.dir.z)
            p_2 = linefit.pos.x**2.0+linefit.pos.y**2.0+linefit.pos.z**2.0
            pd = (linefit.pos.x*linefit.dir.x+linefit.pos.y*linefit.dir.y+linefit.pos.z*linefit.dir.z)
            r_2 = self.vertexRad**2.0
            L = -pd - np.sqrt(pd**2.0-p_2+r_2)
            t_offset = float(frame["TimeShiftedMCPEMap_toffset"].value)
            T0_mc = linefit.time -t_offset + L/self.c
            vertex = dataclasses.I3Position(linefit.pos.x+L*direction.x,linefit.pos.y+L*direction.y,linefit.pos.z+L*direction.z)


        # Set the seed as either true vertex, true dir, or true time.
        if self.group == "dir":
            direction = dataclasses.I3Direction(muon.dir.x,muon.dir.y,muon.dir.z)
        elif self.group == "vert":
            p_2 = muon.pos.x**2.0+muon.pos.y**2.0+muon.pos.z**2.0
            pd = (muon.pos.x*muon.dir.x+muon.pos.y*muon.dir.y+muon.pos.z*muon.dir.z)
            r_2 = self.vertexRad**2.0
            L = -pd - np.sqrt(pd**2.0-p_2+r_2)
            vertex = dataclasses.I3Position(muon.pos.x+L*muon.dir.x,muon.pos.y+L*muon.dir.y,muon.pos.z+L*muon.dir.z)
        elif self.group == "time":
            p_2 = muon.pos.x**2.0+muon.pos.y**2.0+muon.pos.z**2.0
            pd = (muon.pos.x*muon.dir.x+muon.pos.y*muon.dir.y+muon.pos.z*muon.dir.z)
            r_2 = self.vertexRad**2.0
            L = -pd - np.sqrt(pd**2.0-p_2+r_2)
            t_offset = float(frame["TimeShiftedMCPEMap_toffset"].value)
            T0_mc = muon.time -t_offset + L/self.c
        
        VTheta = vertex.theta
        VPhi = vertex.phi

        T0 = GetVertexTime(vertex,direction,data,domsUsed)

        if self.group == "time":
            T0 = T0_mc

        # Minimize using scipy
        def func(x):
            vtheta, vphi, theta, phi, t0 = x
            return qFunctor(vtheta, vphi, theta, phi, t0)
        solution = op.minimize(fun=func, 
                               x0=np.array([VTheta, VPhi, direction.theta, direction.phi, T0]), 
                               method='Nelder-Mead')

        #minimizer = minimizer = Minuit(qFunctor,
         #               t0=T0,
          #              error_t0=100.0,
           #             vtheta=VTheta,
            #            error_vtheta=2.0,
             #           limit_vtheta=(-np.pi,np.pi),
              #          vphi=VPhi,
               #         error_vphi=2.0,
                #        limit_vphi=(-2.0*np.pi,2.0*np.pi),
                 #       phi=direction.phi,
                  #      error_phi=2.0,
                   #     limit_phi=(-2.0*np.pi,2.0*np.pi),
                    #    theta=direction.theta,
                     #   error_theta=2.0,
                      #  limit_theta=(-np.pi,np.pi),
                       # errordef=0.5,
           #             print_level=3
                       # )

        #minimizer.fixed["vtheta"] = True
        #minimizer.fixed["vphi"] = True
        #minimizer.fixed["phi"] = True
        #minimizer.fixed["theta"] = True

        #minimizer.migrad()
        #minimizer.simplex()

        #minimizer.errors["t0"]=100.0
        #minimizer.fixed["vtheta"] = False
        #minimizer.fixed["vphi"] = False
        #minimizer.fixed["phi"] = False
        #minimizer.fixed["theta"] = False

        #minimizer.simplex()
        #minimizer.migrad()

        #solution = minimizer.values
            
        # For likelihood
        #vx = self.vertexRad*np.sin(solution['vtheta'])*np.cos(solution['vphi'])
        #vy = self.vertexRad*np.sin(solution['vtheta'])*np.sin(solution['vphi'])
        #vz = self.vertexRad*np.cos(solution['vtheta'])
        #q = dataclasses.I3Position(vx,vy,vz)
        #phi = solution['phi'] 
        #theta = solution['theta']
        #u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        vx = self.vertexRad*np.sin(solution.x[0])*np.cos(solution.x[1])
        vy = self.vertexRad*np.sin(solution.x[0])*np.sin(solution.x[1])
        vz = self.vertexRad*np.cos(solution.x[0])
        q = dataclasses.I3Position(vx,vy,vz)
        phi = solution.x[3] 
        theta = solution.x[2]
        u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        #print("tau = "+str(solution['tau']))

        #print("post fit")
        #print(str(q.x)+","+str(q.y)+","+str(q.z))

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
        recoParticle.speed = self.c
        recoParticle.pos = q
        recoParticle.time = solution.x[4]

        toc = time.perf_counter() # End timer

        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        #frame[self.output+"_nloglike"] =  dataclasses.I3Double(minimizer.fval)
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(solution.fun)
        frame[self.output+"_time1"] = dataclasses.I3Double(toc - tic)
        frame[self.output+"_seed_llhval"] = dataclasses.I3Double(qFunctor(VTheta, VPhi, direction.theta, direction.phi, T0))
        
        self.PushFrame(frame)    

