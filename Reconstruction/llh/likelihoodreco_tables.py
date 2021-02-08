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
from scipy import special as sp                        # For the Gamma function 
import sys, os
from iminuit import Minuit
import argparse
import math as m

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,vertexrad,pdf,time_lim,dist_lim):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    pmt = []
    time = []
    charge = []
    pdf_tables = pdf

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
    pdf_tables = pdf
    table_time_lim = time_lim
    table_dist_lim = dist_lim
    darkprob = 1e-5

    def GetProbability(distance,time) :

      dist_bin = max(min(int(distance),len(pdf_tables)-1),0)                  
      time_bin = max(min(int(time-table_time_lim[0]),len(pdf_tables[dist_bin])-1),0)
      return pdf_tables[dist_bin][time_bin]+darkprob
    
    # The computations from here on require we find the time and distance of closest approach, d_i,c and t_i,c
    def closestApproach(vtheta, vphi, theta, phi):
        # Compute vec{r} - vec{x}
        vx = vertexRad*np.cos(vphi)*np.sin(vtheta)
        vy = vertexRad*np.sin(vphi)*np.sin(vtheta)
        vz = vertexRad*np.cos(vtheta)
        v =np.array([np.sin(vtheta)*np.cos(vphi),np.sin(vtheta)*np.sin(vphi),np.cos(vtheta)])
        dc = []
        tc = []

        for dom in pulse_series.keys() :                                              
          for pulse in pulse_series[dom] :
            x = geo_doms[dom].position.x - vx
            y = geo_doms[dom].position.y - vy
            z = geo_doms[dom].position.z - vz
            # Compute (\vec{r} - vec{x}) dot \vec{v}
            dotprod = x*v[0] + y*v[1] + z*v[2]
            # Compute the final vector components
            x = x - dotprod*v[0]
            y = y - dotprod*v[1]
            z = z - dotprod*v[2]
            # Compute t_i,c and d_i,c
            dc.append(np.sqrt(x*x + y*y + z*z))
            tc.append(dotprod/c)
        return dc, tc
        
    # Given the linefit and other parameters 
    def computeResiduals(dc, tc, t0):
        t = []
        d = []
        for i in range(len(dc)) :
          # Now we find the time of the photon emission
          _tc = tc[i] - dc[i]/(np.tan(theta_c)*c)
          # The first component of the geometric time
          d.append(dc[i]/np.sin(theta_c)) 
          t_geo = d[-1]/c_n
          # Apply our offset time to find the "true" closest approach time array. Multiply by 1E9 to change to nanoseconds
          _tc += t0
          # The total geometric time
          t_geo = t_geo + _tc 
        # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
          t.append(time[i] - t_geo)
        return d, t

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vtheta, vphi, theta, phi, t0): 
        dc, tc = closestApproach(vtheta, vphi, theta, phi)
        d, t = computeResiduals(dc, tc, t0)

        vx = vertexRad*np.sin(vtheta)*np.cos(vphi)
        vy = vertexRad*np.sin(vtheta)*np.sin(vphi)
        vz = vertexRad*np.cos(vtheta)
        
        prob = []
        for i in range(len(d)) :
          prob.append(GetProbability(d[i],t[i]))

        sum_nloglike = 0.0
        for i in range(len(prob)) :
            sum_nloglike -= charge[i]*np.log(prob[i])
        return sum_nloglike

    return likelihoodFunction

def GetVertexTime(vtheta,vphi,pulse_series,geo_doms,vertexrad):                                 
  # turn PMT locations and time hits into numpy arrays for easier numpy
  # algebra

  mean_t = 0.0                                                            
  mean_x = 0.0                                                            
  mean_y = 0.0                                                            
  mean_z = 0.0                                                            
  sum_charge = 0.0;
                                                                                                                        
  for dom in pulse_series.keys() :                                            
    for pulse in pulse_series[dom] :                                                                          
      mean_t = mean_t + pulse.charge*pulse.time                                 
      mean_x = mean_x + pulse.charge*geo_doms[dom].position.x                                    
      mean_y = mean_y + pulse.charge*geo_doms[dom].position.y                
      mean_z = mean_z + pulse.charge*geo_doms[dom].position.z                
      sum_charge = sum_charge+pulse.charge                                 

  mean_t = mean_t/sum_charge                                          
  mean_x = mean_x/sum_charge                                          
  mean_y = mean_y/sum_charge                                          
  mean_z = mean_z/sum_charge                                          

  vx = vertexrad*np.cos(vphi)*np.sin(vtheta)
  vy = vertexrad*np.sin(vphi)*np.sin(vtheta)                          
  vz = vertexrad*np.cos(vtheta)                                       

  dist_phot = 0.0

  for dom in pulse_series.keys() :                                              
    for pulse in pulse_series[dom] :
      dist_phot += pulse.charge*np.sqrt((mean_x-geo_doms[dom].position.x)**2.0+(mean_y-geo_doms[dom].position.y)**2.0+(mean_z-geo_doms[dom].position.z)**2.0)

  dist_phot = dist_phot/sum_charge                                      
  dist = np.sqrt((vx-mean_x)**2.0+(vy-mean_y)**2.0+(vz-mean_z)**2.0)      
  vertextime = mean_t - dist/0.3- dist_phot/(0.3/1.35)                   
  return vertextime 

class likelihoodreco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("seedtrack","Track to seed fit","linefit")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("vertexRad","Radius to put vertex at",500.)
        self.AddParameter("tablesfile","","")

        self.AddOutBox("OutBox")

    def ReadTables(self) :

      infile = open(self.tablefiles,"r")
      lines = infile.readlines()
      linecount = 0
      self.pdf = [[]]
      xcount = 0
      for line in lines :
        splitline = line.split(",",100)
        if linecount == 0 :
          nx = int(splitline[0].replace("\n",""))
          ny = int(splitline[1].replace("\n",""))
          minx = float(splitline[2].replace("\n",""))
          maxx = float(splitline[3].replace("\n",""))
          miny = float(splitline[4].replace("\n",""))
          maxy = float(splitline[5].replace("\n",""))
        else :
          if xcount == nx :
            self.pdf.append([])
            xcount = 0
          for value in splitline :
            self.pdf[-1].append(float(value.replace("\n","")))
        linecount += 1
      self.time_lim = [miny,maxy]
      self.dist_lim = [minx,maxx]

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.n                                       # light in water
        self.theta_c = np.arccos(1./self.n)                       # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 15.                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data  
        self.tablefiles = self.GetParameter("tablesfile")

        if self.tablefiles == "" :
          self.tablefiles = os.getenv('PONESRCDIR')+"/data/fittertables.dat"

        self.ReadTables()    

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 

        data = frame[self.pulseseries]
        # Clean the data to get rid of repeated events
        #data = clean_data(data)
        linefit = frame[self.seedtrack]
        domsUsed = frame['I3Geometry'].omgeo 

        qFunctor = LikelihoodFunctor(data,domsUsed,self.vertexRad,self.pdf,self.time_lim,self.dist_lim) 
        vr = np.sqrt(linefit.pos.x**2.+linefit.pos.y**2.0+linefit.pos.z**2.0)
        VTheta = np.arccos(linefit.pos.z/vr)
        VPhi = 0.0
        if np.sin(VTheta) != 0.0 :
            #VPhi = np.arccos(linefit.pos.x/(vr*np.sin(VTheta)))
            VPhi = np.arctan2(linefit.pos.y,linefit.pos.x)
        T0 = GetVertexTime(VTheta,VPhi,data,domsUsed,self.vertexRad)
          
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
        recoParticle.speed = self.c
        recoParticle.pos = q
        recoParticle.time = 0
            
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle  
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(minimizer.fval)
        self.PushFrame(frame)    

