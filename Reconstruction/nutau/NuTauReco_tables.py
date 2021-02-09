#!/usr/bin/env python                                                                                                
 
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



def LikelihoodFunctor(data,domsUsed,pdf,time_lim,dist_lim):
                
    pulse_series = data
    geo_doms = domsUsed

    c = 0.299792458                                 # speed of light 
    n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/n                                       # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 557                                       # time parameter that has to be fit using simulations or data 
    pdf_tables = pdf
    table_time_lim = time_lim
    table_dist_lim = dist_lim
    darkprob = 1e-5

    def GetProbability(time,distance) :
  
      if time < table_time_lim[0] or time > table_time_lim[1] :
        return 0.0
      if distance < table_dist_lim[0] or distance > table_dist_lim[1] :
        return 0.0

      dist_bin = max(min(int(distance),len(pdf_tables)-1),0)
      time_bin = max(min(int(time-table_time_lim[0]),len(pdf_tables[dist_bin])-1),0)

      return pdf_tables[time_bin][dist_bin]

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(t0,dt,v1x,v1y,v1z,dtheta,dphi,br):
        v2x = v1x + c*dt*np.cos(dphi)*np.sin(dtheta)
        v2y = v1y + c*dt*np.sin(dphi)*np.sin(dtheta)
        v2z = v1z + c*dt*np.cos(dtheta) 
        darkrate = 1.e-7

        nloglike = 0.0

        for dom in pulse_series.keys() :

            pmt_x = geo_doms[dom].position.x
            pmt_y = geo_doms[dom].position.y
            pmt_z = geo_doms[dom].position.z

            for pulse in pulse_series[dom] :
                t = pulse.time
                charge = pulse.charge
                #vertex 1 
                distance = np.sqrt((v1x-pmt_x)**2.0 + (v1y-pmt_y)**2.0 + (v1z-pmt_z)**2.0) 
                t_trav = distance/c_n
                t_resid = t-t0-t_trav
                prob = br*GetProbability(t_resid,distance)

                #vertex2
                distance = np.sqrt((v2x-pmt_x)**2.0 + (v2y-pmt_y)**2.0 + (v2z-pmt_z)**2.0) 
                t_trav = distance/c_n
                t_resid = t-t0-dt-t_trav
                prob += (1.-br)*GetProbability(t_resid,distance)

                prob += darkrate

                nloglike += -charge*np.log(prob)

        return nloglike

    return likelihoodFunction

class nutaureco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","taufit")
        self.AddParameter("tables","tablesfile","")
        self.AddOutBox("OutBox")

    def ReadTables(self) :

        infile = open(self.tablesfiles,"r")
        lines = infile.readlines()
        linecount = 0
        self.pdf = [[]]
        xcount = 0
        ny = 0
        nx = 0
        minx = 0.0
        maxx = 0.0
        miny = 0.0
        maxy = 0.0
        for line in lines :
            splitline = line.split(",",100)
            if linecount == 0 :
                nx = int(splitline[0].replace("\n",""))
                ny = int(splitline[1].replace("\n",""))
                minx = float(splitline[2].replace("\n",""))
                maxx = float(splitline[3].replace("\n",""))
                miny = float(splitline[4].replace("\n",""))
                maxy = float(splitline[5].replace("\n",""))
                linecount += 1
            else :
                if xcount == ny :
                    self.pdf.append([])
                    xcount = 0
                for value in splitline :
                    self.pdf[-1].append(float(value.replace("\n","")))
                    xcount += 1
            linecount += 1
        self.time_lim = [minx,maxx]
        self.dist_lim = [miny,maxy]

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")

        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.n                             # light in water
        self.theta_c = np.arccos(1./self.n)                  # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 15.                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data   
 
        self.tablesfiles = self.GetParameter("Tables")
        if self.tablesfiles == "" :
          self.tablesfiles = default=os.getenv('PONESRCDIR')+"/data/fittertables.dat"

        self.ReadTables() 

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 
        

        data = frame[self.pulseseries]
        domsUsed = frame['I3Geometry'].omgeo
        sumcharge = 0.0

        _T0  = 0.0
        _V1x = 0.0
        _V1y = 0.0
        _V1z = 0.0

        pulsecount = 0
        
        for dom in data.keys() :

            pmt_x = domsUsed[dom].position.x
            pmt_y = domsUsed[dom].position.y
            pmt_z = domsUsed[dom].position.z

            for pulse in data[dom] :
                    pulsecount += 1
                    _T0 += pulse.time*pulse.charge
                    _V1x += pmt_x*pulse.charge
                    _V1y += pmt_y*pulse.charge
                    _V1z += pmt_z*pulse.charge
                    sumcharge += pulse.charge
        
        if pulsecount < 100 :
          return

        _T0 = _T0/sumcharge
        _V1x = _V1x/sumcharge
        _V1y = _V1y/sumcharge
        _V1z = _V1z/sumcharge

        T0 = _T0
        V1x = _V1x
        V1y = _V1y
        V1z = _V1z

        dT = 0.0
        Dtheta = 0.0
        Dphi = 0.0
        Br = 1.0

        qFunctor = LikelihoodFunctor(data,domsUsed,self.pdf,self.time_lim,self.dist_lim)   
          
        minimizer_single = Minuit(qFunctor, 
                            t0=T0,
                            error_t0=1.0,
                            dt=dT,
                            error_dt=1.0,
                            limit_dt=(0.0,10000.0),
                            v1x=V1x,
                            error_v1x=1.0,
                            limit_v1x=(-500.,500.),
                            v1y=V1y,
                            error_v1y=1.0,
                            limit_v1y=(-500.,500.),
                            v1z=V1z,
                            error_v1z=1.0,
                            limit_v1z=(-500.,500.),
                            dtheta=Dtheta,
                            error_dtheta=1.0,
                            limit_dtheta=(0.0,np.pi),
                            dphi=Dphi,
                            error_dphi=1.0,
                            limit_dphi=(0.0,2.0*np.pi),
                            br=Br,
                            error_br=1.0,
                            limit_br=(0.0,1.0),
                            errordef=0.5,
                           )

        minimizer_single.fixed["dt"]=True
        minimizer_single.fixed["dtheta"]=True
        minimizer_single.fixed["dphi"]=True
        minimizer_single.fixed["br"]=True

        minimizer_single.migrad()

        T0=_T0
        dT = 60*0.3
        V1x = _V1x
        V1y = _V1y
        V1z = _V1z-30
        Dtheta = 0.0
        Dphi = 0.0
        Br = 0.5

        minimizer_double = Minuit(qFunctor,                                     
                                            t0=T0,                                 
                                            error_t0=1.0,
                                            dt=dT,                                              
                                            error_dt=1.0,                                       
                                            limit_dt=(0.0,10000.0),                            
                                            v1x=V1x,                                           
                                            error_v1x=1.0,                                     
                                            limit_v1x=(-500.,500.),                             
                                            v1y=V1y,                                            
                                            error_v1y=1.0,                                      
                                            limit_v1y=(-500.,500.),                           
                                            v1z=V1z,                                           
                                            error_v1z=1.0,                                      
                                            limit_v1z=(-500.,500.),                             
                                            dtheta=Dtheta,                                      
                                            error_dtheta=1.0,                                   
                                            limit_dtheta=(0.0,np.pi),                           
                                            dphi=Dphi,                                          
                                            error_dphi=1.0,                                     
                                            limit_dphi=(0.0,2.0*np.pi),                         
                                            br=Br,                                            
                                            error_br=1.0,                                       
                                            limit_br=(0.0,1.0),                                 
                                            errordef=0.5,                                       
                                            )

        solution_single = minimizer_single.values
        loglikelihood_single = minimizer_single.fval

        v0x = solution_single['v1x']
        v0y = solution_single['v1y']
        v0z = solution_single['v1z']

        v0 = dataclasses.I3Position(v0x,v0y,v0z)
        recoParticle_cascade = dataclasses.I3Particle()
        recoParticle_cascade.shape = dataclasses.I3Particle.InfiniteTrack
        recoParticle_cascade.dir = dataclasses.I3Direction(0.0,0.0,0.0)
        recoParticle_cascade.speed = self.c
        recoParticle_cascade.pos = v0
        recoParticle_cascade.time = solution_single['t0']

        minimizer_double.migrad()

        solution_double = minimizer_double.values
        loglikelihood_double = minimizer_double.fval
            
        # For likelihood
        v1x = solution_double['v1x']
        v1y = solution_double['v1y']
        v1z = solution_double['v1z']

        v2x = v1x + self.c*solution_double['dt']*np.cos(solution_double['dphi'])*np.sin(solution_double['dtheta'])
        v2y = v1y + self.c*solution_double['dt']*np.sin(solution_double['dphi'])*np.sin(solution_double['dtheta'])
        v2z = v1z + self.c*solution_double['dt']*np.cos(solution_double['dtheta'])
      
        v1 = dataclasses.I3Position(v1x,v1y,v1z)
        v2 = dataclasses.I3Position(v2x,v2y,v2z)
        phi = solution_double['dphi'] 
        theta = solution_double['dtheta']
        u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        intensity = solution_double['br']

        # Record the final result
        recoParticle_create = dataclasses.I3Particle()
        recoParticle_create.shape = dataclasses.I3Particle.InfiniteTrack
        recoParticle_create.dir = u
        recoParticle_create.speed = self.c
        recoParticle_create.pos = v1
        recoParticle_create.time = solution_double['t0']

        recoParticle_decay = dataclasses.I3Particle()
        recoParticle_decay.shape = dataclasses.I3Particle.InfiniteTrack                        
        recoParticle_decay.dir = u
        recoParticle_decay.speed = self.c
        recoParticle_decay.pos = v2
        recoParticle_decay.time = solution_double['t0']+solution_double['dt']

            
        # include both linefit and improved recos for comparison
        frame[self.output+"_double_v1"] = recoParticle_create
        frame[self.output+"_double_v2"] = recoParticle_decay 
        frame[self.output+"_double_t0"] = dataclasses.I3Double(intensity)
        frame[self.output+"_double_intensity"] = dataclasses.I3Double(intensity)
        frame[self.output+"_double_nlogl"] =  dataclasses.I3Double(loglikelihood_double)
        frame[self.output+"_single_v0"] = recoParticle_cascade
        frame[self.output+"_single_nlogl"] =  dataclasses.I3Double(loglikelihood_single)
        self.PushFrame(frame)    

