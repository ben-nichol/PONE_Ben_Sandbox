#!/usr/bin/env python                                                                                                
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey 
from icecube.dataclasses import I3Particle 
import numpy as np
import time
from scipy import special as sp                        # For the Gamma function 
from scipy import optimize as op
from scipy import interpolate as interp
from iminuit import Minuit
import sys
#from iminuit import Minuit
import argparse
import math as m
from Utilities.DOMUtility import AddPMTKey
from Utilities.OpticalParameters import c, ngroup, theta_c

# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data,domsUsed,_minpos):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    distance = data
    geo_doms = domsUsed
    minpos = _minpos

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vx,vy,vz):

        vertex = dataclasses.I3Position(vx,vy,vz)
        leastsquare = 0.0
        deltal = np.sqrt((minpos.x-vertex.x)**2.0+(minpos.y-vertex.y)**2.0+(minpos.z-vertex.z)**2.0)
        for dom in distance.keys() :
            dompos = geo_doms[dom].position
            vertex_to_dom  = np.sqrt((dompos.x-vertex.x)**2.0+(dompos.y-vertex.y)**2.0+(dompos.z-vertex.z)**2.0)
            for dist in distance[dom] :
                leastsquare += dist[1]*(vertex_to_dom - dist[0] - deltal)**2.0

        return leastsquare

    return likelihoodFunction

def DirectionFit(data,domsUsed,_vertex):
    distance = data
    geo_doms = domsUsed
    vertex = _vertex
    costheta = np.cos(theta_c)

    def likelihoodFunction(theta,phi):

        direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))

        leastsquare = 0.0
        for dom in distance.keys() :
            dompos = domsUsed[dom].position
            vert_to_dom = dataclasses.I3Direction(dompos.x-vertex.x,dompos.y-vertex.y,dompos.z-vertex.z)
            dot = vert_to_dom.x*direction.x+vert_to_dom.y*direction.y+vert_to_dom.z*direction.z
            for dist in distance[dom] :
                leastsquare += dist[1]*(dot-costheta)**2.0
        return leastsquare

    return likelihoodFunction


def GetDistances(pulse_series,geo_doms,mindoms):                                 

    c_n = c/ngroup                                     # light in water

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
    distance = {}
    for domkey in pulse_series.keys() :
        nopmtkey = AddPMTKey(domkey,1)
        if nopmtkey in mindoms :
            continue
        if nopmtkey not in distance.keys() :
            distance[nopmtkey] = []
        for pulse in pulse_series[domkey] :
            dist = abs(pulse.time - mintime)*c_n
            distance[nopmtkey].append((dist,pulse.charge))

    mindoms.append(mintime_dom)
    
    return T0, minpos, distance

class CascadeReco(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","llnfit")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation. 

    def Geometry(self,frame):

        self.domsUsed = frame['I3Geometry'].omgeo
        self.PushFrame(frame)

    def Physics(self,frame):

        if not frame.Has(self.pulseseries) :
            self.PushFrame(frame)
            return

        data = frame[self.pulseseries]
        
        mindoms = []
        results = []      
        #print("new event")
        for i in range(3) :
            minz = 100000.
            maxz = -100000.
            miny = 100000.
            maxy = -100000.
            minx = 100000.
            maxx = -100000.
            for dom in data.keys() :
                pos = self.domsUsed[dom].position
                minz = min(minz,pos.z)
                maxz = max(maxz,pos.z)
                miny = min(miny,pos.y)
                maxy = max(maxy,pos.y)
                minx = min(minx,pos.x)
                maxx = max(maxx,pos.x)
            minz += -200.
            maxz += 200.
            miny += -200.
            maxy += 200.
            minx += -200.
            maxx += 200. 

            T0, minpos, distance = GetDistances(data,self.domsUsed,mindoms)
            vertex = minpos
            qFunctor = LikelihoodFunctor(distance,self.domsUsed,minpos)

            #print(minpos)
            #print(mindoms)

            # Minimize using scipy
            def func(vx,vy,vz):
                return qFunctor(vx,vy,vz)

            func.errordef = Minuit.LEAST_SQUARES

            m = Minuit(func,
                        vx = vertex.x,
                        vy = vertex.y,
                        vz = vertex.z
                        )

            m.limits["vx"] = (minx,maxx)
            m.limits["vy"] = (miny,maxy)
            m.limits["vz"] = (minz,maxz)
            m.errors["vx"] = 50.0
            m.errors["vy"] = 50.0 
            m.errors["vz"] = 50.0

            m.simplex().migrad()
            res = m.values

            q = dataclasses.I3Position(res["vx"],res["vy"],res["vz"])

            totalcharge = 0.0
            totalcharge_dir = 0.0
            T0 = 0.0
            c_n = c/ngroup
            direc = [0.0,0.0,0.0]

            for dom in data.keys():
                if AddPMTKey(dom,1) in mindoms :
                    continue
                dompos = self.domsUsed[dom].position
                Vertex_to_dom = np.sqrt((dompos.x-q.x)**2.0+(dompos.y-q.y)**2.0+(dompos.z-q.z)**2.0)
                if Vertex_to_dom > 0.0 :
                    _direct = [(dompos.x-q.x)/Vertex_to_dom,(dompos.y-q.y)/Vertex_to_dom,(dompos.z-q.z)/Vertex_to_dom]
                for pulse in data[dom] :
                    direc[0] += pulse.charge*_direct[0]
                    direc[1] += pulse.charge*_direct[1]
                    direc[2] += pulse.charge*_direct[2]
                    T0 += pulse.charge*(pulse.time - Vertex_to_dom/c_n)
                    totalcharge += pulse.charge
                    if Vertex_to_dom > 0.0 :
                        totalcharge_dir += pulse.charge
            T0 /= totalcharge
            direc[0] /= totalcharge_dir
            direc[1] /= totalcharge_dir
            direc[2] /= totalcharge_dir

            seeddir  = dataclasses.I3Direction(direc[0],direc[1],direc[2])

            qFunctor2 = DirectionFit(distance,self.domsUsed,q)

            def func2(theta,phi):
                return qFunctor2(theta,phi)

            func2.errordef = Minuit.LEAST_SQUARES

            m2 = Minuit(func2,
                        theta = seeddir.theta,
                        phi = seeddir.phi
                        )

            m2.limits["theta"] = (0.0,np.pi)
            m2.limits["phi"] = (0.0,2.0*np.pi)
            m2.errors["theta"] = 1.0
            m2.errors["phi"] = 1.0

            m2.simplex().migrad()
            res2 = m2.values

            direction = dataclasses.I3Direction(np.sin(res2["theta"]*np.cos(res2["phi"])),np.sin(res2["theta"])*np.sin(res2["phi"]),np.cos(res2["theta"]))

            # Record the final result
            recoParticle = dataclasses.I3Particle()
            recoParticle.shape = dataclasses.I3Particle.Cascade

            # record on particle whether reconstruction was successful
            #if solution.success == True:
            #    recoParticle.fit_status = dataclasses.I3Particle.OK
            #else:
            #    recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality
                                            
            recoParticle.speed = c
            recoParticle.pos = q
            recoParticle.dir = direction
            recoParticle.time = T0 

            results.append((recoParticle,m.fval/len(distance)))
            # include both linefit and improved recos for comparison

        finalresult = results[0][0]
        finalval = results[0][1]
        finalfit = 0
        for i in range(len(results)) :
            #print(results[0][1])
            #print(results[0][0])
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

