#!/usr/bin/env python                                                                                                
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey 
from icecube.dataclasses import I3Particle 
import numpy as np
from iminuit import Minuit
import sys
#from iminuit import Minuit
import argparse
import math as m
from Utilities.DOMUtility import AddPMTKey, DOMProperties
from Utilities.OpticalParameters import c,GetIndex,GetGroupIndex ,theta_c


# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def MinimizerFunctor(data,domsUsed,_minpos):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    deltaT = data
    geo_doms = domsUsed
    minpos = _minpos
    c_ngroup = c/GetGroupIndex()
    sintheta = np.sin(theta_c)
    tantheta = np.tan(theta_c)

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def minimizerFunction(vx,vy,vz,theta,phi):

        direct = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))

        vertex = dataclasses.I3Position(vx,vy,vz)
        
        T0 = -np.sqrt((minpos.x-vertex.x)**2.0+(minpos.y-vertex.y)**2.0+(minpos.z-vertex.z)**2.0)/c_ngroup

        leastsquare = 0.0
        for dom in deltaT.keys() :
            dompos = geo_doms[dom].position
            vertex_to_dom  = dataclasses.I3Position(dompos.x-vertex.x,dompos.y-vertex.y,dompos.z-vertex.z) 
            dc = vertex_to_dom.x*direct.x+vertex_to_dom.y*direct.y+vertex_to_dom.z*direct.z
            l_c = dataclasses.I3Position(vertex_to_dom.x - dc*direct.x,vertex_to_dom.y - dc*direct.y,vertex_to_dom.z - dc*direct.z)
            l_photon = l_c.r/sintheta
            l_emmission = dc-l_c.r/tantheta

            np.sqrt((dompos.x-vertex.x)**2.0+(dompos.y-vertex.y)**2.0+(dompos.z-vertex.z)**2.0)
            for delt in deltaT[dom] :
                leastsquare += delt[1]*(delt[0]-T0 - l_emmission/c - l_photon/c_ngroup)**2.0

        return leastsquare

    return minimizerFunction

def GetDeltaT(pulse_series,geo_doms,mindoms):                                 

    mintime = 10000.
    mintime_dom = None
    
    for domkey in pulse_series.keys() :
        if AddPMTKey(domkey,1) in mindoms :
            continue
        for pulse in pulse_series[domkey] :
            if pulse.time < mintime :
                mintime = pulse.time
                mintime_dom = AddPMTKey(domkey,1)

    T0 = mintime

    minpos = geo_doms[mintime_dom].position
    deltaT = {}
    Photondir_x = 0.0
    Photondir_y = 0.0
    Photondir_z = 0.0
    totalcharge = 0.0

    for domkey in pulse_series.keys() :
        nopmtkey = AddPMTKey(domkey,1)
        if nopmtkey in mindoms :
            continue
        if nopmtkey not in deltaT.keys() :
            deltaT[nopmtkey] = []
        for pulse in pulse_series[domkey] :
            pulse.time - mintime
            deltaT[nopmtkey].append((pulse.time - mintime,pulse.charge))
            totalcharge += pulse.charge
            
    mindoms.append(mintime_dom)
    
    return T0, minpos, deltaT


class TrackReco_leastsquare(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddOutBox("OutBox")
        self.AddParameter("vertexRad","The Radius to readjust the vertex to",200.)
        self.AddParameter("seedtrack","name of seed track","linefit")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.seedtrack = self.GetParameter("seedtrack")

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo

        maxradius = 0.0
        for dom in self.domsUsed.keys() :
            pos = self.domsUsed[dom].position
            radius = np.sqrt(pos.x**2.0+pos.y**2.0+pos.z**2.0)
            maxradius = max(maxradius,radius)

        self.vertexRad += maxradius

        self.PushFrame(frame)

    def Physics(self,frame):

        if not frame.Has(self.pulseseries) :
            self.PushFrame(frame)
            return

        data = frame[self.pulseseries]
        if self.seedtrack != None :
            linefit = frame[self.seedtrack].dir
        else :
            linefit = dataclasses.I3Direction(0.0,0.0,1.0)

        mindoms = []
        results = []
        for i in range(3) :
            T0, minpos, deltaT = GetDeltaT(data,self.domsUsed,mindoms)
            vertex = minpos
            qFunctor = MinimizerFunctor(deltaT,self.domsUsed,minpos)

            def func(vx,vy,vz,theta,phi):
                return qFunctor(vx,vy,vz,theta,phi)

            func.errordef = Minuit.LEAST_SQUARES

            m = Minuit(func,
                        vx = minpos.x,
                        vy = minpos.y,
                        vz = minpos.z,
                        theta = linefit.theta,
                        phi = linefit.phi,
                        )

            m.limits["vx"] = (minpos.x-200.,minpos.x+200.)
            m.limits["vy"] = (minpos.y-200.,minpos.y+200.)
            m.limits["vz"] = (minpos.z-200.,minpos.z+200.)
            m.limits["theta"] = (0.0,np.pi)
            m.limits["phi"] = (0.0,2.0*np.pi)
            m.errors["vx"] = 50.0
            m.errors["vy"] = 50.0
            m.errors["vz"] = 50.0
            m.errors["theta"] = 1.0
            m.errors["phi"] = 1.0

            m.simplex().migrad()

            res = m.values

            c_ngroup = c/GetGroupIndex()

            dx = np.sin(res["theta"])*np.cos(res["phi"])
            dy = np.sin(res["theta"])*np.sin(res["phi"])
            dz = np.cos(res["theta"])

            direct = dataclasses.I3Direction(dx,dy,dz)

            vx = res["vx"]
            vy = res["vy"]
            vz = res["vz"]

            T0 += -np.sqrt((minpos.x -vx)**2.0+(minpos.y -vy)**2.0+(minpos.z -vz)**2.0)/c_ngroup

            dot = vx*dx + vy*dy + vz*dz
            vv = vx*vx + vy*vy + vz*vz
            rr = self.vertexRad*self.vertexRad
            
            sqr2 = dot*dot+rr-vv
            l = 0.0
            if sqr2 > 0.0 :
                l = -dot - np.sqrt(sqr2)
            elif sqr2 == 0.0 :
                l = -dot

            vx += l*dx
            vy += l*dy
            vz += l*dz

            T0 += -l/c
            
            q = dataclasses.I3Position(vx,vy,vz)
            direction = dataclasses.I3Direction(dx,dy,dz)

            recoParticle = dataclasses.I3Particle()
            recoParticle.shape = dataclasses.I3Particle.Cascade
                                            
            recoParticle.speed = c
            recoParticle.pos = q
            recoParticle.dir = direction
            recoParticle.time = T0 

            results.append((recoParticle,m.fval/len(deltaT)))

        finalresult = results[0][0]
        finalval = results[0][1]
        finalfit = 0
        for i in range(len(results)) :
            if results[i][1] < finalval :
                finalresult = results[i][0]
                finalval = results[i][1]
                finalfit = i

        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        for dom in data.keys() :
            nopmtkey = AddPMTKey(dom,1)
            if finalfit > 0 and nopmtkey ==  mindoms[0] :
                continue
            if finalfit > 1 and nopmtkey == mindoms[1] :
                continue
            outputpulsemap[dom] = data[dom]

        frame[self.output] = finalresult
        frame[self.output+"_leastsquare"] =  dataclasses.I3Double(finalval)
        frame[self.output+"_pulseseries"] = outputpulsemap
        self.PushFrame(frame)    

