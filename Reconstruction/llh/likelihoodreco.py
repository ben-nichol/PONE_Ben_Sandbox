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
from Reconstruction.llh.reco_pdfs import cpandel as pdf               # This module is used to store the pdf
from scipy import special as sp                        # For the Gamma function 
import sys
from iminuit import Minuit
import argparse
import math as m

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,vertexrad):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    pmt = []
    time = []
    charge = []

    if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
        for dom in pulse_series.keys() :
            for pulse in pulse_series[dom] :
                pmt.append([])
                pmt[-1].append(geo_doms[dom].position.x)
                pmt[-1].append(geo_doms[dom].position.y)
                pmt[-1].append(geo_doms[dom].position.z)
                time.append(pulse.time)
                charge.append(pulse.charge)
    else :
        for dom in pulse_series.keys() :
            for pulse in pulse_series[dom] :
                pmt.append([])
                pmt[-1].append(geo_doms[dom].position.x)
                pmt[-1].append(geo_doms[dom].position.y)
                pmt[-1].append(geo_doms[dom].position.z)
                time.append(pulse.time)
                charge.append(1.0)  

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
    min_index = np.argmin(time)
    
    # The computations from here on require we find the time and distance of closest approach, d_i,c and t_i,c
    def closestApproach():
        # Compute vec{r} - vec{x}
        dc = []
        tc = []
        lc = []

        for dom in pulse_series.keys() :                                              
          for pulse in pulse_series[dom] :
            x = geo_doms[dom].position.x - vx
            y = geo_doms[dom].position.y - vy
            z = geo_doms[dom].position.z - vz
            # Compute (\vec{r} - vec{x}) dot \vec{v}
            dotprod = x*v[0] + y*v[1] + z*v[2]
            # Compute the final vector components
            # Compute t_i,c and d_i,c
            dc.append(np.sqrt(x*x + y*y + z*z-dotprod*dotprod))
            #dc.append(np.sqrt(x*x + y*y + z*z))
            tc.append(dotprod/c)
            lc.append(dotprod)
        return dc, tc, lc
        
    # Given the linefit and other parameters 
    def computeResiduals(dc, tc, lc, t0):
        t = []
        d = []
        l = []
        for i in range(len(dc)) :
          # Now we find the time of the photon emission
          _tc = tc[i] - dc[i]/(np.tan(theta_c)*c)
          l.append(lc[i]-dc[i]/np.tan(theta_c))
          # The first component of the geometric time
          d.append(dc[i]/np.sin(theta_c)) 
          t_geo = d[-1]/c_n
          # Apply our offset time to find the "true" closest approach time array. Multiply by 1E9 to change to nanoseconds
          _tc += t0
          # The total geometric time
          t_geo = t_geo + _tc 
        # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
          t.append(time[i] - t_geo)
        return d, t, l

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vtheta, vphi, theta, phi, t0):
        vx = vertexRad*np.sin(vtheta)*np.cos(vphi)
        vy = vertexRad*np.sin(vtheta)*np.sin(vphi)
        vz = vertexRad*np.cos(vtheta)
        v =np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        dc, tc, lc = closestApproach()
        d, t, l = computeResiduals(dc, tc, lc, t0)
        p_charge=[]
        dark = 1.e-8

        N = 0.0
        for dom in pulse_series.keys():
          totalcharge = 0.0
          for pulse in pulse_series[dom]:
             totalcharge += pulse.charge
          x = geo_doms[dom].position.x - vx
          y = geo_doms[dom].position.y - vy
          z = geo_doms[dom].position.z - vz
          # Compute (\vec{r} - vec{x}) dot \vec{v}
          dotprod = x*v[0] + y*v[1] + z*v[2]
          # Compute the final vector components
          # Compute t_i,c and d_i,c
          d_c = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
          l = dotprod - d_c/np.tan(theta_c)
          d_p = d_c/np.sin(theta_c)
          N += totalcharge/(np.exp(-d_p/tau)/max(d_c,0.25)+1000.*dark)


        for i in range(len(dc)) :
        #  if l[i]<start or l[i]>start+stop :
        #    p_charge.append(0.0)
        #  else :
        #    p_charge.append(np.exp(-d[i]/tau)/max(dc[i],0.5))
            p_charge.append(np.exp(-d[i]/tau)/max(dc[i],0.25))
        out = pdf(t,d)
        sum_nloglike = 0.0
        for i in range(len(out)) :
            sum_nloglike -= charge[i]*np.log(out[i]*p_charge[i]+dark/N)

	#likelihood of not seeing any light
	#assume chargetotal = sum(N*p_charge[i]) -> N = chargetotal/sum(p_charge[i])
        #Hit=[]
        #NotHit=[]
        #for dom in geo_doms.keys() :





          #x = geo_doms[dom].position.x - vx
          #y = geo_doms[dom].position.y - vy
          #z = geo_doms[dom].position.z - vz
          # Compute (\vec{r} - vec{x}) dot \vec{v}
          #dotprod = x*v[0] + y*v[1] + z*v[2]
          # Compute the final vector components
          #x = x - dotprod*v[0]
          #y = y - dotprod*v[1]
          #z = z - dotprod*v[2]
          # Compute t_i,c and d_i,c
          #d_c = np.sqrt(x*x + y*y + z*z)
          #l = dotprod - d_c/np.tan(theta_c) 
          #d_p = d_c/np.sin(theta_c)
          #if l < start or l > start+stop :
          #  if dom in pulse_series.keys() :
          #    Hit.append(1000.*dark)
          #  else :
          #    NotHit.append(1000.*dark)
          #else :
          #  if dom in pulse_series.keys() :
          #    Hit.append(np.exp(-d_p/tau)/max(d_c,0.5)+1000.*dark)
          #  else :
          #    NotHit.append(np.exp(-d_p/tau)/max(d_c,0.5)+1000.*dark)

        #N = sum(charge)/sum(Hit)
        #for p in NotHit :
        #  sum_nloglike += p*N 

        return sum_nloglike

    return likelihoodFunction

def GetStartStop(vertex,direction,pulse_series,geo_doms) :
        n = 1.34
        theta_c = np.arccos(1./n)  

        start = 99999999999
        stop = 0.0

        for dom in pulse_series.keys() :
          dompos = geo_doms[dom].position

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
  # turn PMT locations and time hits into numpy arrays for easier numpy
  # algebra

	totalcharge = 0.0
	MaxChargeDOM = 0
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
				MaxChargeDOM = dom
	else :
		for dom in pulse_series.keys() :
			for pulse in pulse_series[dom] :
				totalcharge += 1.0
			if totalcharge > maxCharge :
                                maxCharge = totalcharge
                                MaxChargeDOM = dom

	#time of largest pulse
	maxCharge=0.0
	maxCharge_time = 0.0
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
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.useMC = self.GetParameter("UseMC")

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

        data = frame[self.pulseseries]
        # Clean the data to get rid of repeated events
        #data = clean_data(data)
	
        linefit = frame[self.seedtrack]
        domsUsed = frame['I3Geometry'].omgeo

        direction = dataclasses.I3Direction(linefit.dir.x,linefit.dir.y,linefit.dir.z)
        
        #linefit = frame[self.seedtrack]
        domsUsed = frame['I3Geometry'].omgeo 

        qFunctor = LikelihoodFunctor(data,domsUsed,self.vertexRad)

        p_2 = linefit.pos.x**2.0+linefit.pos.y**2.0+linefit.pos.z**2.0
        pd = (linefit.pos.x*direction.x+linefit.pos.y*direction.y+linefit.pos.z*direction.z)
        r_2 = self.vertexRad**2.0

        if pd**2.0-p_2+r_2 < 0.0 :
            return

        L = -pd - np.sqrt(pd**2.0-p_2+r_2)

        vertex = dataclasses.I3Position(linefit.pos.x+L*direction.x,linefit.pos.y+L*direction.y,linefit.pos.z+L*direction.z)

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
            
        VTheta = vertex.theta
        VPhi = vertex.phi

        T0 = GetVertexTime(vertex,direction,data,domsUsed)

        minimizer = minimizer = Minuit(qFunctor,
                        t0=T0,
                        error_t0=100.0,
                        vtheta=VTheta,
                        error_vtheta=2.0,
                        limit_vtheta=(-np.pi,np.pi),
                        vphi=VPhi,
                        error_vphi=2.0,
                        limit_vphi=(-2.0*np.pi,2.0*np.pi),
                        phi=direction.phi,
                        error_phi=2.0,
                        limit_phi=(-2.0*np.pi,2.0*np.pi),
                        theta=direction.theta,
                        error_theta=2.0,
                        limit_theta=(-np.pi,np.pi),
                        errordef=0.5,
           #             print_level=3
                        )

        minimizer.fixed["vtheta"] = True
        minimizer.fixed["vphi"] = True
        minimizer.fixed["phi"] = True
        minimizer.fixed["theta"] = True

        minimizer.migrad()

        minimizer.errors["t0"]=100.0
        minimizer.fixed["vtheta"] = False
        minimizer.fixed["vphi"] = False
        minimizer.fixed["phi"] = False
        minimizer.fixed["theta"] = False

        minimizer.migrad()

        solution = minimizer.values
            
        # For likelihood
        vx = self.vertexRad*np.sin(solution['vtheta'])*np.cos(solution['vphi'])
        vy = self.vertexRad*np.sin(solution['vtheta'])*np.sin(solution['vphi'])
        vz = self.vertexRad*np.cos(solution['vtheta'])
        q = dataclasses.I3Position(vx,vy,vz)
        phi = solution['phi'] 
        theta = solution['theta']
        u = dataclasses.I3Direction(-np.sin(theta)*np.cos(phi), -np.sin(theta)*np.sin(phi), -np.cos(theta))

        #print("tau = "+str(solution['tau']))

        #print("post fit")
        #print(str(q.x)+","+str(q.y)+","+str(q.z))

        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.InfiniteTrack
                
        # record on particle whether reconstruction was successful
        if minimizer.get_fmin()["is_valid"]:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
        recoParticle.dir = u
        recoParticle.speed = self.c
        recoParticle.pos = q
        recoParticle.time = solution['t0']
            
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(minimizer.fval)
        self.PushFrame(frame)    

