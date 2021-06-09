#!/usr/bin/env python                                                                                                
# This is meant to be a slightly more robust approach to reconstruction of a muon event.                               
# The physics and likelihood model is heavily based off of the ICECUBE model and can be found at                      
# "https://publications.ub.uni-mainz.de/theses/volltexte/2014/3869/pdf/3869.pdf"                                     
# The time residuals are computed by myself though. The techniques used are detailed in a text document I have somewhere -dg
 
# Import some useful ICECUBE modules                                                                                  
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey  
from icecube.dataclasses import I3Particle 
from Reconstruction.llh.reco_pdfs import cpandel as pdf               # This module is used to store the pdf
from scipy import special as sp                        # For the Gamma function 
import sys
from iminuit import Minuit
from scipy.optimize import minimize, LinearConstraint, basinhopping
import argparse
import math as m
import numpy as np

def GetPMTAcceptanceTable(infile) :

    domaccFile = open(infile,"r")
    lines = domaccFile.readlines()
    PMTacceptance = list()

    zenithcount = 0
    for line in lines :

        splitline = line.split(" ",1000)
        if zenithcount % 179 == 0 :
            zenithcount = 0
            PMTacceptance.append([])
        PMTacceptance[-1].append([])
        for value in splitline :
            PMTacceptance[-1][-1].append(float(value))
        zenithcount += 1
    return PMTacceptance

def GetPMTAcceptance(PMTacceptance,photonDir,pmtid):

    if pmtid > len(PMTacceptance)-1 :
         print("Bad PMTID")
         print(pmtid)

    thetaBin = max(0,min(178,int(180.0*photonDir.zenith/np.pi)))
    phiBin = max(0,min(358,int(180.0*photonDir.azimuth/np.pi)))
    if thetaBin > len(PMTacceptance[pmtid]) -1 :
        print("Bad thetaBin")
        print(thetaBin)
    if phiBin > len(PMTacceptance[pmtid][thetaBin]) -1 :
        print("Bad thetaBin")
        print(phiBin)
    return PMTacceptance[pmtid][thetaBin][phiBin]

def ComputeCentroid(pulse_series,geo_doms) :

    centroid_x = 0.0
    centroid_y = 0.0
    centroid_z = 0.0
    allcharge = 0.0

    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        pos = geo_doms[domkey].position

        totalcharge = 0.0
        for pulse in pulse_series[dom] :
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                totalcharge += pulse.charge
            else :
                totalcharge  += 1.0
        
        centroid_x += pos.x*totalcharge
        centroid_y += pos.y*totalcharge
        centroid_z += pos.z*totalcharge
        allcharge += totalcharge

    centroid_x /= allcharge
    centroid_y /= allcharge
    centroid_z /= allcharge

    return dataclasses.I3Position(centroid_x,centroid_y,centroid_z)

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def GetGeoTime(position,vert) :
        c = 0.299792458                                 # speed of light 
        n = 1.34
        ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
        c_n = c/ngroup                                     # light in water
        x = position.x - vert.x
        y = position.y - vert.y
        z = position.z - vert.z
        d = m.sqrt(x*x + y*y + z*z)
        t = d/c_n 
        return d,t

def sumloglike(t0,vertex,pulse_series,geo_doms,PMTacceptance):
    dark = 1.e-8
    tau = 18.949132224466762

    sum_nloglike = 0.0
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, dom.pmt)
        d,t = GetGeoTime(geo_doms[domkey].position,vertex)
        p_charge = np.exp(-d/tau)/max(d**2.0,0.25**2.0)
        pmtaccept = 1.0 #GetPMTAcceptance(PMTacceptance,pdir,dom.pmt)
        for pulse in pulse_series[dom] :
            charge = 1.0
            cpandel_out = pdf(pulse.time - t0 - t ,d)
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            sum_nloglike -= charge*np.log(cpandel_out*p_charge*pmtaccept+dark)
            sum_nloglike -= charge*min(0.0,pulse.time - t0 - t)
    return sum_nloglike

def TimeLikelihoodFunctor(data,domsUsed,_vertex,_PMTacceptance) :
 # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    tau = 18.949132224466762
    vertex = _vertex
    PMTacceptance = _PMTacceptance
    def likelihoodFunction(t0):
        return sumloglike(t0,vertex,pulse_series,geo_doms,PMTacceptance)

    return likelihoodFunction

def LikelihoodFunctor(data,domsUsed,vrad,_PMTacceptance):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    vertex_radius = vrad
    PMTacceptance = _PMTacceptance
    T0 = 6000

    def likelihoodFunction(rad,theta,z):
    #def likelihoodFunction(x):
        vx = rad*np.cos(theta)
        vy = rad*np.sin(theta)
        vz = z
        vertex = dataclasses.I3Position(vx,vy,vz)
        timeFunctor = TimeLikelihoodFunctor(pulse_series,geo_doms,vertex,PMTacceptance)
        minimizer = Minuit(timeFunctor,
                              t0=T0,
                              error_t0=100.0,
                              errordef=0.5,
                              #print_level=4
                             )
        minimizer.migrad()
        #T0 = minimizer.values['t0']
        #func = sumloglike(T0,vertex,direction,pulse_series,geo_doms,PMTacceptance)
        return minimizer.fval

    return likelihoodFunction

class nutaureco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("seedtrack","Track to seed fit","linefit")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddParameter("vertexRad","Radius to put vertex at",550.)
        self.AddParameter("UseMC","Use MC Truth Track to seed",False)
        self.AddParameter("SplitDoms","",False)
        self.AddParameter("DOMAcceptanceFile","","")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.useMC = self.GetParameter("UseMC")
        self.splitDOMs = self.GetParameter("SplitDoms")
        self.mDOMPMTaccept = list()
        if self.GetParameter("DOMAcceptanceFile") != "" :
             self.mDOMPMTaccept = GetPMTAcceptanceTable(self.GetParameter("DOMAcceptanceFile"))
        else :
             self.splitDOMs = False

        # Some quantities that are environment dependent
        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34  
        self.ngroup = 1.3555714017                                      # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.ngroup                                       # light in water
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 18.949132224466762                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data      

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 

        data = frame[self.pulseseries]
	
        direction = dataclasses.I3Direction(0.0,0.0,0.0)
        
        domsUsed = frame['I3Geometry'].omgeo 

        T0 = 6000

        vertex = ComputeCentroid(data,domsUsed) 

        prefitParticle = dataclasses.I3Particle()
        prefitParticle.shape = dataclasses.I3Particle.InfiniteTrack
        prefitParticle.dir = direction
        prefitParticle.speed = self.c
        prefitParticle.pos = vertex
        prefitParticle.time = T0

        finalFunctor = LikelihoodFunctor(data,domsUsed,self.vertexRad,self.mDOMPMTaccept)
        minimizer_final = Minuit(finalFunctor,
                        rad = vertex.rho,
                        error_rad = 1.0,
                        limit_rad = (1.0,200.0),
                        theta = vertex.phi,
                        error_theta = 0.5,
                        limit_theta = (vertex.phi-0.5,vertex.phi+0.5),
                        z = vertex.z,
                        error_z = 1.0,
                        limit_z = (-500.0,500.0),
                        errordef = 0.5,
                        )

        minimizer_final.migrad()
        solution_final = minimizer_final.values
        likelihood = minimizer_final.fval

        solution_initial = solution_final	

        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.Cascade
                                    
        vx = solution_final['rad']*np.cos(solution_final['theta'])
        vy = solution_final['rad']*np.sin(solution_final['theta'])
        vz = solution_final['z']
        vertex = dataclasses.I3Position(vx,vy,vz)

        timeFunctor = TimeLikelihoodFunctor(data,domsUsed,vertex,self.mDOMPMTaccept)
        minimizer_t = Minuit(timeFunctor,
                              t0=6000.,
                              error_t0=100.0,
                              errordef=0.5,
                             )
        minimizer_t.migrad()
        T0 = minimizer_t.values['t0']
        
        recoParticle.dir = direction
        recoParticle.speed = self.c
        recoParticle.pos = vertex
        recoParticle.time = T0
    
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle 
        frame[self.output+"_prefit"] = prefitParticle 
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(likelihood)
        self.PushFrame(frame)    

