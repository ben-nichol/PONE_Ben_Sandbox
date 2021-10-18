#!/usr/bin/env python

'''

'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from numpy import linalg as la
from icecube.phys_services import I3Calculator
import operator
from Utilities.DOMUtility import NoPMTKey, AddPMTKey

class CausalPulseCleaning(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("inputseries","GCD to be simulated","MCPESeriesMap")
        self.AddParameter("output","GCD to be simulated","causalpulses")
        self.AddParameter("window","Window used to decide most sig window.",200)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.gcdFile = self.GetParameter("GCDFile")
        self.input = self.GetParameter("inputseries")
        self.output = self.GetParameter("output")
        self.window = self.GetParameter("window")

    #Compute base time and position that pulses need to be causal 
    #with if no other pulses are causal with it.
    def getBaseTimeandPosition(self,mcpeList):
        DOMCharge = {}

        for omkey in mcpeMap.keys():
            DOMkey = NoPMTKey(omkey)
            if DOMkey not in DOMCharge.keys():
                DOMCharge[DOMkey] = 0.0
            for pulse in mcpeMap[omkey] :
                DOMCharge[DOMkey] += pulse.charge
            
        pulses = list()
        MaxchargeDOM = DOMCharge.keys()[0]
        for DOM in DOMCharge.keys():
            if DOMCharge[DOM] > DOMCharge[MaxchargeDOM] :
                MaxchargeDOM = DOM

            npmts = GetNPMTs()
            for ipmt in range(npmts):
                PMTKey = OMKey(MaxchargeDOM.string,MaxchargeDOM.om,ipmt)
                if PMTkey in mcpeMap.keys():
                    for pulse in mcpeMap[PMTKey] :
                        pulses.append(pulses.time)
        pulsecoinc = list()
        for ipulse1 in range(len(pulses)) :
                pulsecoinc.append(0)
                for ipulse2 in range(len(pulses)) :
                    if abs(pulses[ipulse1]-pulses[ipulse2]) < 10. :
                        pulsecoinc[-1] += 1

        maxcoince = max(pulsecoinc)
        coincetime = max(pulses)

        for i in range(len(pulsecoince)):
            if pulsecoinc[i] == maxcoince :
                if pulses[i] < coincetime :
                    coincetime = pulses[i]

        return MaxchargeDOM, coincetime

    def getCausalMCPEs(self,pos1,pos2,mcpeList,coincetime):
        causalMCPEList = dataclasses.I3RecoPulseSeries()

        dis = np.sqrt((pos1.x-pos2.x)**2.0+(pos1.y-pos2.y)**2.0+(pos1.z-pos2.z)**2.0)

        windowmin = coincetime - dis*1.5/0.3
        windowmax = coincetime + dis*1.5/0.3

        for pulse in mcpeList:
          if pulse.time >  windowmin and pulse.time <  windowmax :
            causalMCPEList.append(pulse)

        return causalMCPEList
   
    def DAQ(self,frame):
    
        mcpeMap = frame[self.input]
        causalMCPEMap = dataclasses.I3RecoPulseSeriesMap()
        domsUsed = frame['I3Geometry'].omgeo

        #make assumption that DOM with highest charge is the baseline, 
        #the 10 ns window in theis DOM with the highest charge is 
        # the baseline for other photons being causal.

        MaxchargeDOM, coincetime = self.getBaseTimeandPosition(self,mcpeMap)

        pos1 = domsUsed[MaxChargeDOM].position
        for omkey in mcpeMap.keys():
            #if omkey not in mcpeMap:
            #  continue
            if len(mcpeMap[omkey]) < 1 :
              continue
            pos2 = domsUsed[NoPMTKey(omkey)].position
            causalHist = self.getCausalMCPEs(pos1,pos2,mcpeMap[omkey],coincetime)
            significantMCPEMap[omkey] = self.getSignificantMCPEs(mcpeMap[omkey])
    
        frame[self.output] = causalMCPEMap
        self.PushFrame(frame)
