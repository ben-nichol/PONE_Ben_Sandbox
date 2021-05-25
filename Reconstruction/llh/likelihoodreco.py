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

def GetVertex(vertexRad,intvertex,direction) :

    p_2 = intvertex.x**2.0+intvertex.y**2.0+intvertex.z**2.0
    pd = (intvertex.x*direction.x+intvertex.y*direction.y+intvertex.z*direction.z)
    r_2 = vertexRad**2.0

    if pd**2.0-p_2+r_2 < 0.0 :
        return dataclasses.I3Position(0.0,0.0,vertexRad)

    L1 = -pd -  np.sqrt(pd**2.0-p_2+r_2)
    L2 = -pd +  np.sqrt(pd**2.0-p_2+r_2)
    v1 = dataclasses.I3Position(intvertex.x+L1*direction.x,intvertex.y+L1*direction.y,intvertex.z+L1*direction.z)
    v2 = dataclasses.I3Position(intvertex.x+L2*direction.x,intvertex.y+L2*direction.y,intvertex.z+L2*direction.z)

    return v1,v2  

def GetVertexPoint(direction,centroid) :
   center = centroid
   direct = direction
   def computePoint(p1,p2) :
      if direct.x> max(direct.y,direct.z) :
         p3 = (direct.y*p1+direct.z*p2)/direct.x
         return dataclasses.I3Position(p3+center.x,p1+center.y,p2+center.z)
      if direct.y > direct.z :
         p3 = (direct.x*p1+direct.z*p2)/direct.y
         return dataclasses.I3Position(p1+center.x,p3+center.y,p2+center.z)
      p3 = (direct.x*p1+direct.y*p2)/direct.z
      return dataclasses.I3Position(p1+center.x,p2+center.y,p3+center.z)

   return computePoint

def ComputeCentroid(vertex,direction,pulse_series,geo_doms) :

     chargesegment = np.zeros(2000)
     meantime = np.zeros(2000)
     counts = np.zeros(2000)
     c = 0.299792458                                 # speed of light 
     ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
     c_n = c/ngroup  

     for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        t -= d/c_n
        #print(t) 
        binval  = max(0,min(int(t*c),1999))
        #print(binval) 
        for pulse in pulse_series[dom] :
            meantime[binval] += pulse.time - d/c_n
            counts[binval] += 1.0
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 chargesegment[binval] += pulse.charge
            else :
                 chargesegment[binval]  += 1.0

     for i in range(len(meantime)) :
         if counts[i]>0.0 :
             meantime[i] = meantime[i]/counts[i]
     lengths = []
     time = []
     maxcharge = 0.0
     L = 0.0
     for i in range(len(chargesegment)) :
         L += i*chargesegment[i]
     L = L/sum(chargesegment)
     return dataclasses.I3Position(vertex.x+L*direction.x,vertex.y+L*direction.y,vertex.z+L*direction.z)


# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
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
        dc = 0.25
        if (x*x + y*y + z*z-dotprod*dotprod) > 0.5 :
            dc = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
        d = dc/np.sin(theta_c)
        perp = dc/np.tan(theta_c)-dotprod
        pdir = dataclasses.I3Direction(direction.x*perp+position.x,direction.y*perp+position.y,direction.z*perp+position.z) 
        t = d/c_n - perp/c 
        return d,dc,t,pdir

def sumloglike(t0,vertex,direction,pulse_series,geo_doms,PMTacceptance):
    dark = 1.e-8
    tau = 18.949132224466762

    sum_nloglike = 0.0
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, dom.pmt)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        p_charge = np.exp(-d/tau)/max(dc,0.25)
        pmtaccept = 1.0 #GetPMTAcceptance(PMTacceptance,pdir,dom.pmt)
        first = True
        for pulse in pulse_series[dom] :
            charge = 1.0
            cpandel_out = pdf(pulse.time - t0 - t ,d)
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            sum_nloglike -= charge*np.log(cpandel_out*p_charge*pmtaccept+dark)
            sum_nloglike -= charge*min(0.0,pulse.time - t0 - t)
    return sum_nloglike

def chargesumloglike(vertex,direction,pulse_series,geo_doms,PMTacceptance):
    dark = 1.e-8
    tau = 18.949132224466762

    sum_nloglike = 0.0
    N = 0.0
    ndoms = 0.0
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        p_charge = np.exp(-d/tau)/max(dc,0.25)
        pmtaccept = 1.0 #GetPMTAcceptance(PMTacceptance,pdir,dom.pmt)
        totalcharge = 0.0
        ndoms += 1.0
        for pulse in pulse_series[dom] :
            charge = 1.0
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            totalcharge += charge
        N+= totalcharge / p_charge*pmtaccept
        #sum_nloglike -= totalcharge*np.log(p_charge+dark)
    N = N/ndoms

    for dom in geo_doms.keys() :
        for i in range(len(PMTacceptance)):
            domkey =  OMKey(dom.string, dom.om, i)
            if domkey in pulse_series.keys():
                continue
            d,dc,t,pdir = GetGeoTime(geo_doms[dom].position,vertex,direction)
            p_charge = np.exp(-d/tau)/max(dc,0.25)
            pmtaccept = 1.0 #GetPMTAcceptance(PMTacceptance,pdir,domkey.pmt)
            sum_nloglike += N*p_charge*pmtaccept

    return sum_nloglike


def timeresidualvariance(vertex,direction,pulse_series,geo_doms):
    dark = 1.e-8
    tau = 18.949132224466762
    sigma = 1.0
    d_o = 50.0
    a_o = 100.00
    d_1 = 0.25

    sum_nloglike = 0.0
    npulses = 0.0
    ndoms = 0.0
    sumcharge = 0.0
    sum_timeres = 0.0
    sum_chargefrac = 0.0
    chargefracs = []
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        totalcharge = 0.0
        for pulse in pulse_series[dom] :
            charge = 1.0
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            totalcharge += charge
            npulses += 1.0
        A = totalcharge*a_o/np.sqrt(a_o**2.0+totalcharge**2.0)
        D = np.sqrt(d_1**2.0+dc**2.0)
        sum_chargefrac += A*D
        chargefracs.append(A*D)
        sumcharge += A
        ndoms += 1.0

    meanchargefrac = sum_chargefrac/ndoms
    sigmachargefrac = 0.0
    for i in range(len(chargefracs)) :
        sigmachargefrac += (chargefracs[i]-meanchargefrac)**2.0

    return sum_chargefrac/((sumcharge/ndoms)*d_o)

def timesumloglike(t0,vertex,direction,pulse_series,geo_doms):
    dark = 1.e-8
    tau = 18.949132224466762

    sum_nloglike = 0.0
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        p_charge = np.exp(-d/tau)/max(dc,0.25)
        for pulse in pulse_series[dom] :
            charge = 1.0
            cpandel_out = pdf(pulse.time - t0 - t ,d)
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            sum_nloglike -= charge*np.log(cpandel_out+dark)
            sum_nloglike -= charge*min(0.0,pulse.time - t0 - t)
    return sum_nloglike

def TimeLikelihoodFunctor(data,domsUsed,_vertex,_direction,_PMTacceptance) :
 # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    tau = 18.949132224466762
    vertex = _vertex
    direction = _direction
    PMTacceptance = _PMTacceptance
    def likelihoodFunction(t0):
        return sumloglike(t0,vertex,direction,pulse_series,geo_doms,PMTacceptance)

    return likelihoodFunction

def LikelihoodFunctor(data,domsUsed,vrad,_PMTacceptance):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    vertex_radius = vrad
    PMTacceptance = _PMTacceptance
    T0 = 6000

    def likelihoodFunction(vtheta,vphi, theta, phi):
    #def likelihoodFunction(x):
        vx = vertex_radius*np.sin(vtheta)*np.cos(vphi)
        vy = vertex_radius*np.sin(vtheta)*np.sin(vphi)
        vz = vertex_radius*np.cos(vtheta)
        vertex = dataclasses.I3Position(vx,vy,vz)
        dx = np.sin(theta)*np.cos(phi)
        dy = np.sin(theta)*np.sin(phi)
        dz = np.cos(theta)
        direction = dataclasses.I3Direction(dx,dy,dz)
        timeFunctor = TimeLikelihoodFunctor(pulse_series,geo_doms,vertex,direction,PMTacceptance)
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

def VarianceFunctor(data,domsUsed,vradius):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    tau = 18.949132224466762
    vertex_r = vradius
 
    def likelihoodFunction(vtheta,vphi,theta,phi):
    #def likelihoodFunction(x):
        vx = vertex_r*np.sin(vtheta)*np.cos(vphi)
        vy = vertex_r*np.sin(vtheta)*np.sin(vphi)
        vz = vertex_r*np.cos(vtheta)
        vertex = dataclasses.I3Position(vx,vy,vz)
        dx = np.sin(theta)*np.cos(phi)
        dy = np.sin(theta)*np.sin(phi)
        dz = np.cos(theta)
        direction = dataclasses.I3Direction(dx,dy,dz)
        func = timeresidualvariance(vertex,direction,pulse_series,geo_doms)
        return func

    return likelihoodFunction

def ChargeLikelihoodFunctor(data,domsUsed,vertexrad,_PMTacceptance):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    tau = 18.949132224466762
    PMTacceptance = _PMTacceptance
    vertexRad = vertexrad

    def likelihoodFunction(vtheta, vphi, theta, phi):
    #def likelihoodFunction(x):
        vx = vertexRad*np.sin(vtheta)*np.cos(vphi)
        vy = vertexRad*np.sin(vtheta)*np.sin(vphi)
        vz = vertexRad*np.cos(vtheta)
        vertex = dataclasses.I3Position(vx,vy,vz)
        dx = np.sin(theta)*np.cos(phi)
        dy = np.sin(theta)*np.sin(phi)
        dz = np.cos(theta)
        direction = dataclasses.I3Direction(dx,dy,dz)

        return chargesumloglike(vertex,direction,pulse_series,geo_doms,PMTacceptance)

    return likelihoodFunction

class likelihoodreco(icetray.I3ConditionalModule):

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

        vertex,vertex2 = GetVertex(self.vertexRad,linefit.pos,direction)
        T0 = 6000

        #centroid,T0 = ComputeCentroid(vertex,direction,data,domsUsed) 

        if self.useMC :
        #if True :
            MMCTrackList = frame['MMCTrackList']
            muon = MMCTrackList[0].GetI3Particle()
            direction_muon = dataclasses.I3Direction(muon.dir.x,muon.dir.y,muon.dir.z)
            p_2 = muon.pos.x**2.0+muon.pos.y**2.0+muon.pos.z**2.0
            pd = (muon.pos.x*muon.dir.x+muon.pos.y*muon.dir.y+muon.pos.z*muon.dir.z)
            r_2 = self.vertexRad**2.0
            L = -pd - np.sqrt(pd**2.0-p_2+r_2)
            t_offset = float(frame["TimeShiftedMCPEMap_toffset"].value)
            T0_mc = linefit.time -t_offset + L/self.c
            vx = muon.pos.x+L*direction_muon.x
            vy = muon.pos.y+L*direction_muon.y
            vz = muon.pos.z+L*direction_muon.z
            vertex_muon = dataclasses.I3Position(vx,vy,vz)
        #    print("truth")
        #    print(vertex_muon)
        #    print(direction_muon)

        #print("linefit")
        #print(vertex)
        #print(direction) 

        #vFunctor = VarianceFunctor(data,domsUsed,self.vertexRad)
        #qFunctor = ChargeLikelihoodFunctor(data,domsUsed,self.vertexRad,self.mDOMPMTaccept)
        #minimizer = Minuit(qFunctor,
        #                vtheta = vertex.theta,
        #                error_vtheta = 0.25,
        #                limit_vtheta = (vertex.theta-0.25,vertex.theta+0.25),
        #                vphi = vertex.phi,
        #                error_vphi = 0.25,
        #                limit_vphi = (vertex.phi-0.25,vertex.phi+0.25),
        #                phi = direction.phi,
        #                error_phi = 0.25,
        #                limit_phi = (direction.phi-0.25,direction.phi+0.25),
        #                theta = direction.theta,
        #                error_theta = 0.25,
        #                limit_theta = (direction.theta-0.25,direction.theta+0.25),
        #                errordef = 0.5,
                        #print_level=4
        #                )
        #print(minimizer.values)
        #minimizer.migrad()
        #solution = minimizer.values
        #print("solution")
        #print(solution)

        #likelihood = minimizer.fval
        #solution_final = solution
   
        #vx = self.vertexRad*np.sin(solution_final['vtheta'])*np.cos(solution_final['vphi'])
        #vy = self.vertexRad*np.sin(solution_final['vtheta'])*np.sin(solution_final['vphi'])
        #vz = self.vertexRad*np.cos(solution_final['vtheta'])
        #vertex = dataclasses.I3Position(vx,vy,vz)

        #dx = np.sin(solution_final['theta'])*np.cos(solution_final['phi'])
        #dy = np.sin(solution_final['theta'])*np.sin(solution_final['phi'])
        #dz = np.cos(solution_final['theta'])
        #direction = dataclasses.I3Direction(dx,dy,dz)
        #print("first fit")
        #print(direction)
        #print(vertex)


        #print("likelihood")
        #print(likelihood)
        #centroid = ComputeCentroid(vertex,direction,data,domsUsed)
        #print("centroid")
        #print(centroid)
        #for i in range(10) :
        #    for j in range(20) :
        #        _theta = np.pi/2.0+((linefit.dir.theta-np.pi/2.0)/np.abs(linefit.dir.theta-np.pi/2.0))*float(i)*(np.pi/20.)
        #        _phi = float(j)*(np.pi/10.)
        #        dx = np.sin(_theta)*np.cos(_phi)
        #        dy = np.sin(_theta)*np.sin(_phi)
        #        dz = np.cos(_theta)
        #        direction = dataclasses.I3Direction(dx,dy,dz)
        #        vertex,vertex2 = GetVertex(self.vertexRad,centroid,direction)
        #        minimizer_final = Minuit(qFunctor,
        #                                         vtheta = vertex.theta,
        #                                         error_vtheta = 0.25,
        #                                         limit_vtheta = (0.0,np.pi),
        #                                         vphi = vertex.phi,
        #                                         error_vphi = 0.25,
        #                                         limit_vphi = (0.0,2.0*np.pi),
       #                                          phi = direction.phi,
       #                                          error_phi = 0.25,
       #                                          limit_phi = (0.0,2.0*np.pi),
       #                                          theta = direction.theta,
       #                                          error_theta = 0.25,
       #                                          limit_theta = (0.0,np.pi),
       #                                          errordef = 0.5,
       #                          #               print_level=4
       #                                         )
       #         minimizer_final.migrad()
       #         if minimizer_final.fval < likelihood :
       #             solution_final = minimizer_final.values
       #             likelihood = minimizer_final.fval
       # print("grid search")
       # print(solution_final)
       # print("likelihood")
       # print(likelihood)

        #vertexfunc = GetVertexPoint(direction,centroid)
        #vertexinner = vertexfunc(solution['p1'],solution['p2'])
        #vx = self.vertexRad*np.sin(solution_final['vtheta'])*np.cos(solution_final['vphi'])
        #vy = self.vertexRad*np.sin(solution_final['vtheta'])*np.sin(solution_final['vphi'])
        #vz = self.vertexRad*np.cos(solution_final['vtheta'])
        #vertex = dataclasses.I3Position(vx,vy,vz)

        #dx = np.sin(solution_final['theta'])*np.cos(solution_final['phi'])
        #dy = np.sin(solution_final['theta'])*np.sin(solution_final['phi'])
        #dz = np.cos(solution_final['theta'])
        #direction = dataclasses.I3Direction(dx,dy,dz)
        #print(direction)
        #print(vertex)
        #vertex,vertex2 = GetVertex(self.vertexRad,vertex,direction)
        #print(vertex)
        #timeFunctor = TimeLikelihoodFunctor(data,domsUsed,vertex,vertex2,self.mDOMPMTaccept)
        #minimizer_t1 = Minuit(timeFunctor,
        #                      t0=6000,
        #                      error_t0=100.0,
        #                      errordef=0.5,
                              #print_level=4
        #                     )
        #minimizer_t1.migrad()

        #timeFunctor = TimeLikelihoodFunctor(data,domsUsed,vertex2,vertex,self.mDOMPMTaccept)
        #minimizer_t2 = Minuit(timeFunctor,
        #                      t0=6000,
        #                      error_t0=100.0,
        #                      errordef=0.5,
			      #print_level=4
        #                     )
        #minimizer_t2.migrad()
        
        
        #T0 = 6000
        #if minimizer_t1.fval < minimizer_t2.fval :
        #T0 = minimizer_t1.values['t0']
            #plrint("use T0_1")
        #else :
        #    T0 = minimizer_t2.values['t0']
        #    direction = dataclasses.I3Direction(vertex.x-vertex2.x,vertex.y-vertex2.y,vertex.z-vertex2.z)
        #    vertex = vertex2
            #print("use T0_2")
        #print(vertex)
        prefitParticle = dataclasses.I3Particle()
        prefitParticle.shape = dataclasses.I3Particle.InfiniteTrack
        prefitParticle.dir = direction
        prefitParticle.speed = self.c
        prefitParticle.pos = vertex
        prefitParticle.time = T0

        finalFunctor = LikelihoodFunctor(data,domsUsed,self.vertexRad,self.mDOMPMTaccept)
        minimizer_final = Minuit(finalFunctor,
         #               t0 = T0,
         #               error_t0 = 5.0,
                        vtheta = vertex.theta,
                        error_vtheta = 0.5,
                        limit_vtheta = (vertex.theta-0.5,vertex.theta+0.5),
                        vphi = vertex.phi,
                        error_vphi = 0.5,
                        limit_vphi = (vertex.phi-0.5,vertex.phi+0.5),
                        phi = direction.phi,
                        error_phi = 0.5,
                        limit_phi = (direction.phi-0.5,direction.phi+0.5),
                        theta = direction.theta,
                        error_theta = 0.5,
                        limit_theta = (direction.theta-0.5,direction.theta+0.5),
                        errordef = 0.5,
#         #               print_level=4
                        )
        #print(minimizer_final.values)
        minimizer_final.migrad()
        solution_final = minimizer_final.values
        likelihood = minimizer_final.fval
        #print("solutions")
        #print(solution_final)
        solution_initial = solution_final	

	#find the charge centroid on the line and pick vertexes and directions to intersect.
#        centroid = ComputeCentroid(vertex,direction,data,domsUsed) 
#        for i in range(200) :	
#            T0 = solution_initial['t0']-1000.+10.*i
#            minimizer_final = Minuit(finalFunctor,
#                                                 t0 = T0,
#                                                 error_t0 = 10.0,
#                                                 vtheta = solution_initial['vtheta'],
#                                                 error_vtheta = 0.25,
#                                                 limit_vtheta = (solution_initial['vtheta']-0.25,solution_initial['vtheta']+0.25),
#                                                 vphi = solution_initial['vphi'],
#                                                 error_vphi = 0.25,
#                                                 limit_vphi = (solution_initial['vphi']-0.25,solution_initial['vphi']+0.25#),
#                                                 phi = solution_initial['phi'],
#                                                 error_phi = 0.25,
#                                                 limit_phi = (solution_initial['phi']-0.25,solution_initial['phi']+0.25),
#                                                 theta = solution_initial['theta'],
#                                                 error_theta = 0.25,
#                                                 limit_theta = (solution_initial['theta']-0.25,solution_initial['theta']+0.25),
#                                                 errordef = 0.5,
                                #               print_level=4
#                                                )
#            minimizer_final.fixed['t0'] = True
#            minimizer_final.migrad()
#            minimizer_final.fixed['t0'] = False
#            minimizer_final.migrad()
#            if m.isnan(minimizer_final.fval) or m.isnan(minimizer_final.values['vtheta']) :
#                continue
            #print("t0")
            #print(T0)
            #print(minimizer_final.fval)
            #print(minimizer_final.values)
#            if minimizer_final.fval < likelihood :
#                solution_final = minimizer_final.values
#                likelihood = minimizer_final.fval
#
#        print("solutions")
#        print(solution_final) 
        #T0 = solution_final['t0'] 	

        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.InfiniteTrack
                                    
        vx = self.vertexRad*np.sin(solution_final['vtheta'])*np.cos(solution_final['vphi'])
        vy = self.vertexRad*np.sin(solution_final['vtheta'])*np.sin(solution_final['vphi'])
        vz = self.vertexRad*np.cos(solution_final['vtheta'])
        vertex = dataclasses.I3Position(vx,vy,vz)

        dx = np.sin(solution_final['theta'])*np.cos(solution_final['phi'])
        dy = np.sin(solution_final['theta'])*np.sin(solution_final['phi'])
        dz = np.cos(solution_final['theta'])
        direction = dataclasses.I3Direction(dx,dy,dz)

        timeFunctor = TimeLikelihoodFunctor(data,domsUsed,vertex,direction,self.mDOMPMTaccept)
        minimizer_t = Minuit(timeFunctor,
                              t0=6000.,
                              error_t0=100.0,
                              errordef=0.5,
                              #print_level=4
                             )
        minimizer_t.migrad()
        T0 = minimizer_t.values['t0']

        #print(vertex)
        #print(direction)
        #T0 = solution_final['t0']	
        
        recoParticle.dir = direction
        recoParticle.speed = self.c
        recoParticle.pos = vertex
        recoParticle.time = T0

    
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle 
        frame[self.output+"_prefit"] = prefitParticle 
        frame[self.output+"_nloglike"] =  dataclasses.I3Double(likelihood)
        self.PushFrame(frame)    

