#!/usr/bin/env python

'''

'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import numpy as np
from numpy import linalg as la
from icecube.phys_services import I3Calculator
import operator
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, GetNPMTs
from Utilities.OpticalParameters import c, n, ngroup, theta_c

class CausalPulseCleaning(icetray.I3ConditionalModule):
    '''
    Causal Pulse Cleaning looks at the brightest DOM and assigns a time and location such that all other 
    hits must be within a loose causality window. This Window is based on the travel of the particle at 
    c + at most 50m of photon travel time at velocity c/n.
    '''
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("inputseries","GCD to be simulated","MCPESeriesMap")
        self.AddParameter("output","GCD to be simulated","causalpulses")
        self.AddParameter("windowscale","Scale the causality window by factor to make it a looser cut",1.2)
        self.AddParameter("maxphotonatten","Maximum number of attenuationlengths photons can travel",2.3)
        self.AddParameter("MaxTracktoReferenceDOM","Maximum distance between Reference DOM and track",50.0)
        self.AddParameter("attenlength","Atteneuation length of the light",35.0)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.gcdFile = self.GetParameter("GCDFile")
        self.input = self.GetParameter("inputseries")
        self.output = self.GetParameter("output")
        self.windowscale = self.GetParameter("windowscale")
        self.maxphotonatten = self.GetParameter("maxphotonatten")
        self.MaxTracktoReferenceDOM = self.GetParameter("MaxTracktoReferenceDOM")
        self.attenlength = self.GetParameter("attenlength")
    #Compute base time and position that pulses need to be causal 
    #with if no other pulses are causal with it.
    def getBaseTimeandPosition(self,mcpeMap):
        DOMCharge = {}

        for omkey in mcpeMap.keys():
            DOMkey = NoPMTKey(omkey)
            if DOMkey not in DOMCharge.keys():
                DOMCharge[DOMkey] = 0.0
            for pulse in mcpeMap[omkey] :
                DOMCharge[DOMkey] += pulse.charge
            
        pulses = list()
        first = True
        MaxchargeDOM = None 
        for DOM in DOMCharge.keys():
            if first :
                MaxchargeDOM = DOM
                first = False
            if DOMCharge[DOM] > DOMCharge[MaxchargeDOM] :
                MaxchargeDOM = DOM

        npmts = GetNPMTs()
        for ipmt in range(npmts):
            PMTkey = AddPMTKey(MaxchargeDOM,ipmt)
            if PMTkey in mcpeMap.keys():
                for pulse in mcpeMap[PMTkey] :
                    pulses.append(pulse.time)

        pulsecoinc = list()
        for ipulse1 in range(len(pulses)) :
                pulsecoinc.append(0)
                for ipulse2 in range(len(pulses)) :
                    if abs(pulses[ipulse1]-pulses[ipulse2]) < 10. :
                        pulsecoinc[-1] += 1

        maxcoince = max(pulsecoinc)
        coincetime = max(pulses)

        for i in range(len(pulsecoinc)):
            if pulsecoinc[i] == maxcoince :
                if pulses[i] < coincetime :
                    coincetime = pulses[i]

        return MaxchargeDOM, coincetime

    def getCausalMCPEs(self,pos1,pos2,mcpeList,coincetime):
        causalMCPEList = dataclasses.I3RecoPulseSeries()

        dis = np.sqrt((pos1.x-pos2.x)**2.0+(pos1.y-pos2.y)**2.0+(pos1.z-pos2.z)**2.0)


        dis_maxphoton = min(dis,self.attenlength*self.maxphotonatten)
        dis_part = 0.0
        if dis_maxphoton < dis :
            dis_part = np.sin(theta_c-np.arcsin(np.sin(np.pi-theta_c)*(self.attenlength*self.maxphotonatten/dis)))*dis/np.sin(np.pi-theta_c)

        WindowMin = max(0.0,(dis/c)-self.MaxTracktoReferenceDOM*ngroup/c)/self.windowscale
        WindowMax = self.windowscale*(dis_part/c + dis_maxphoton*ngroup/c)

        if pos1==pos2 :
            WindowMin = -300.
            WindowMax = 300.

        for pulse in mcpeList:
            deltaT = abs(pulse.time-coincetime)
            if deltaT > WindowMin and deltaT < WindowMax :
                causalMCPEList.append(pulse)

        return causalMCPEList
   
    def DAQ(self,frame):
    
        mcpeMap = frame[self.input]
        causalMCPEMap = dataclasses.I3RecoPulseSeriesMap()
        domsUsed = frame['I3Geometry'].omgeo

        #make assumption that DOM with highest charge is the baseline, 
        #the 10 ns window in theis DOM with the highest charge is 
        # the baseline for other photons being causal.

        MaxchargeDOM, coincetime = self.getBaseTimeandPosition(mcpeMap)

        pos1 = domsUsed[MaxchargeDOM].position
        for omkey in mcpeMap.keys():
            #if omkey not in mcpeMap:
            #  continue
            if len(mcpeMap[omkey]) < 1 :
              continue
            pos2 = domsUsed[NoPMTKey(omkey)].position
            causalHits = self.getCausalMCPEs(pos1,pos2,mcpeMap[omkey],coincetime)
            if len(causalHits) > 0 :
                causalMCPEMap[omkey] = causalHits
    
        frame[self.output] = causalMCPEMap
        frame[self.output+"_string"] = dataclasses.I3Double(MaxchargeDOM.string)
        frame[self.output+"_om"] = dataclasses.I3Double(MaxchargeDOM.om)
        frame[self.output+"_coincetime"] = dataclasses.I3Double(coincetime)
        self.PushFrame(frame)
