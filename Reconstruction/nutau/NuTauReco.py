#!/usr/bin/env python                                                                                                
 
# Import some useful ICECUBE modules                                                                                  
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
from scipy import special as sp                        # For the Gamma function 
import sys
from iminuit import Minuit
import argparse
import math as m

def amplitude(dist) :

    dist = max(dist,1.0)
    absorbtion = 50.0
    lambda_s = 120.
    atten = np.sqrt(lambda_s**2.0 + absorbtion**2.0)
    return np.exp(-dist/atten)*(1./dist**2.0)

def cpandel(t, d, sigma = 2.0, lambda_s = 120., rho = 0.004):
    xi = d/lambda_s
    eta = rho*sigma - (t/sigma)   

    if xi <= 5. and (rho*sigma - (300)/sigma <= eta and eta < (rho*sigma+5.) ) :
    # Define our region dependent approximations of the CPandel function
        first = xi*np.log(rho) + (xi - 1.)*np.log(sigma) - (t**2)/(2.*sigma**2) - ((1. + xi)/2.)*np.log(2)
        frac_1 = sp.hyp1f1(0.5*xi,0.5,0.5*eta**2)/sp.gamma(0.5*(xi + 1.))
        frac_2 = sp.hyp1f1(0.5*(xi+1.),1.5,0.5*eta**2)/sp.gamma(0.5*xi)
        second = np.log(frac_1 - np.sqrt(2.)*eta*frac_2)
        return -(first + second)

    if xi <= 1. and eta < (rho*sigma - (300)/sigma) :
        first = (rho**2)*(sigma**2)/2.
        second = xi*np.log(rho) + (xi - 1.)*np.log(t) - np.log(xi) - rho*t - sp.gamma(xi)
        return -(first + second)  

    if (xi>5. and eta<0.) or (xi>1. and eta<(rho*sigma-(300)/sigma)) :
        z = -eta/np.sqrt(4*xi - 2.)
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
        phi = 1. - N1/(2.*xi - 1.) + N2/((2.*xi - 1.)**2)
        alpha = -t**2/(2*sigma**2) + 0.25*eta**2 - xi*0.5 + 0.25 + k*(2*xi - 1.) - 0.25*np.log(1 + z**2) - 0.5*xi*np.log(2) + 0.5*(xi - 1.)*np.log(2*xi - 1.) + xi*np.log(rho) + (xi - 1.)*np.log(sigma)
        return -(alpha - np.log(sp.gamma(xi)) + np.log(phi))  

    if (xi>5. and  eta>= 0.) or (xi>1. and eta>=(rho*sigma+5.)) :
        z = eta/np.sqrt(4*xi - 2.)
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
        psi = 1. + N1/(2*xi - 1.) + N2/((2*xi - 1.)**2)
        first = xi*np.log(rho) + (xi - 1.)*np.log(sigma) - (t**2)/(2*sigma**2) + 0.25*(eta**2) - 0.25*np.log(2*np.pi)
        U = 0.5*xi - 0.25 - 0.5*xi*np.log(2*xi - 1.) + 0.5*(xi - 1.)*np.log(2)
        second = -k*(2*xi - 1.) - 0.25*np.log(1. + z**2) + np.log(psi)
        return -(first + U + second)

    if xi<=1. and eta>=(rho*sigma + 5.) :
        return -(xi*np.log(rho*sigma) - 0.5*np.log(2*np.pi*sigma**2) -
            xi*np.log(eta) - (t**2)/(2*sigma**2)) 

    print("error, no condition met")
    return 0.0

class nutaureco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("GCDFile","GCD file.")
        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","taufit")

        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")
        self.gcdfile = self.GetParameter("GCDFile")
        self.geometry = self.gcdfile.pop_frame()["I3Geometry"]
        self.domsUsed = self.geometry.omgeo

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.n                                       # light in water
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 15.                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data      

    # Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
    def LikelihoodFunctor(self,data, prnt = False):

                
        pulse_series = data

        # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
        def likelihoodFunction(t0,dt,v1x,v1y,v1z,dtheta,dphi,br):
            v2x = v1x + self.c*dt*np.cos(dphi)*np.sin(dtheta)
            v2y = v1y + self.c*dt*np.sin(dphi)*np.sin(dtheta)
            v2z = v1z + self.c*dt*np.cos(dtheta) 
            darkrate = 1./10000.

            nloglike = 0.0

            for dom in pulse_series.keys() :

                pmt_x = self.domsUsed[dom].position.x
                pmt_y = self.domsUsed[dom].position.y
                pmt_z = self.domsUsed[dom].position.z

                for pulse in pulse_series[dom] :
                    t = pulse.time
                    charge = pulse.charge
                    #vertex 1 
                    distance = np.sqrt((v1x-pmt_x)**2.0 + (v1y-pmt_y)**2.0 + (v1z-pmt_z)**2.0) 
                    amp =  amplitude(distance)
                    t_trav = distance/(self.c_n)
                    prob = br*amp*cpandel(t-to-t_trav,distance)

                    #vertex2
                    distance = np.sqrt((v2x-pmt_x)**2.0 + (v2y-pmt_y)**2.0 + (v2z-pmt_z)**2.0) 
                    amp =  amplitude(distance)
                    t_trav = distance/(self.c_n)
                    prob += (1.-br)*amp*cpandel(t-to-dt-t_trav,distance)

                    nloglike += -charge*np.log(prob+darkrate)

            return nloglike

            return likelihoodFunction

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 
        

        data = frame[self.pulseseries]

        sumcharge = 0.0

        T0 = 0.0
        V1x = 0.0
        V1y = 0.0
        V1z = 0.0

        for dom in data.keys() :
            pmt_x = self.domsUsed[dom].position.x
            pmt_y = self.domsUsed[dom].position.y
            pmt_z = self.domsUsed[dom].position.z

            for pulse in data[dom] :
                    T0 += pulse.time*pulse.charge
                    V1x += pmt_x*pulse.charge
                    V1y += pmt_y*pulse.charge
                    V1z += pmt_z*pulse.charge
                    sumcharge += pulse.charge
        
        T0 = T0/sumcharge
        V1x = V1x/sumcharge
        V1y = V1y/sumcharge
        V1z = V1z/sumcharge

        dT = 0.0
        Dtheta = 0.0
        Dphi = 0.0
        Br = 1.0

        qFunctor = self.LikelihoodFunctor(data)   
          
        minimizer = Minuit(qFunctor, 
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
                            print_level=3
                           )

        minimizer.fixed["dt"]=True
        minimizer.fixed["dtheta"]=True
        minimizer.fixed["dphi"]=True
        minimizer.fixed["br"]=True

        minimizer.migrad()

        solution_single = minimizer.values
        loglikelihood_single = minimizer.fval

        v0x = solution_single['v1x']
        v0y = solution_single['v1y']
        v0z = solution_single['v1z']

        v0 = dataclasses.I3Position(v0x,v0y,v0z)
        recoParticle_cascade = dataclasses.I3Particle()
        recoParticle_cascade.shape = dataclasses.I3Particle.InfiniteTrack
        recoParticle_cascade.dir = dataclasses.I3Direction(0.0,0.0,0.0)
        recoParticle_cascade.speed = c
        recoParticle_cascade.pos = v0
        recoParticle_cascade.time = solution_double['t0']

        minimizer.fixed["dt"]=False
        minimizer.fixed["dtheta"]=False
        minimizer.fixed["dphi"]=False
        minimizer.fixed["br"]=False

        minimizer.migrad()

        solution_double = minimizer.values
        loglikelihood_double = minimizer.fval
            
        # For likelihood
        v1x = solution_double['v1x']
        v1y = solution_double['v1y']
        v1z = solution_double['v1z']

        v2x = v1x + self.c*dt*np.cos(solution_double['dphi'])*np.sin(solution_double['dtheta'])
        v2y = v1y + self.c*dt*np.sin(solution_double['dphi'])*np.sin(solution_double['dtheta'])
        v2z = v1z + self.c*dt*np.cos(solution_double['dtheta'])
      
        v1 = dataclasses.I3Position(v1x,v1y,v1z)
        v2 = dataclasses.I3Position(v2x,v2y,v2z)
        phi = solution_double['dphi'] 
        theta = solution_double['dtheta']
        u = dataclasses.I3Direction(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))

        intensity = solutions_double['br']

        # Record the final result
        recoParticle_create = dataclasses.I3Particle()
        recoParticle_create.shape = dataclasses.I3Particle.InfiniteTrack
        recoParticle_create.dir = u
        recoParticle_create.speed = c
        recoParticle_create.pos = v1
        recoParticle_create.time = solution_double['t0']

        recoParticle_decay = dataclasses.I3Particle()
        recoParticle_decay.shape = dataclasses.I3Particle.InfiniteTrack                        
        recoParticle_decay.dir = u
        recoParticle_decay.speed = c
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
        self.push_frame(frame)    

