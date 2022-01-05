from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from math import sqrt
from copy import deepcopy

class DetectorTrigger(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("output","Append the outputs",'')
        self.AddParameter("input","Name of the Physics I3MCTree name","")
        self.AddParameter("CutOnTrigger","Cut events that do not trigger.",False)
        self.AddParameter("FullDetectorCoincidenceN","",3)
        self.AddParameter("FullDetectorCoincidenceWindow","",1.1)
        self.AddParameter("StringCoincidenceN","",2)
        self.AddParameter("StringCoincidenceWindow","",1.1)
        self.AddParameter("StringNRows","",3)
        self.AddParameter("StringDist","",1.5)
        self.AddParameter("ScaleBySpacing","Turn the Windows inputs to be relative to detector size.",True)
        self.AddParameter("ForceAdjacency","Require adjacency ",True)
        self.AddParameter("DOMPMTCoinc","Number of PMTs needed for DOM Trigger",2)
        self.AddParameter("EventLength","Length of Event",10000)
        self.AddParameter("TriggerTime","Time of trigger in event.",2000)
        self.AddParameter("PulseSeriesIn","Pulse series in","")
        self.AddParameter("PulseSeriesOut","Pulse series out","")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.gcdFile = self.GetParameter("GCDFile")
        self.output = self.GetParameter("output")
        self.input = self.GetParameter("input")
        self.CutOnTrigger = self.GetParameter("CutOnTrigger")
        self.FullDetectorCoincidenceN = self.GetParameter("FullDetectorCoincidenceN")
        self.StringCoincidenceN = self.GetParameter("StringCoincidenceN")
        self.ForceAdjacency = self.GetParameter("ForceAdjacency")
        self.StringNRows = self.GetParameter("StringNRows")
        self.StringDist = self.GetParameter("StringDist")
        self.FullDetectorCoincidenceWindow = self.GetParameter("FullDetectorCoincidenceWindow")
        self.StringCoincidenceWindow = self.GetParameter("StringCoincidenceWindow")
        self.DOMPMTCoinc = self.GetParameter("DOMPMTCoinc")
        self.EventLength = self.GetParameter("EventLength")
        self.TriggerTime = self.GetParameter("TriggerTime")
        self.PulseSeriesIn = self.GetParameter("PulseSeriesIn")
        self.PulseSeriesOut = self.GetParameter("PulseSeriesOut")
        self.DoStringTrigger = (self.StringCoincidenceN < self.FullDetectorCoincidenceN) and (self.StringCoincidenceWindow < self.FullDetectorCoincidenceWindow)

        self.nstrings = int(0)
        self.nDOMs = int(0)
        self.npmts = int(0)

        if self.GetParameter("ScaleBySpacing") :
            #Figure out largest distance between two DOMs, average distance between closest String, Average min distance between DOMs.

            #there is a better way for this.
            for frame in self.gcdFile:
                domsUsed = frame['I3Geometry'].omgeo
            self.gcdFile.rewind()
            for domkey in domsUsed.keys() :
                if domkey.string > self.nstrings :
                    self.nstrings = int(domkey.string)
                if domkey.om > self.nDOMs  :
                    self.nDOMs  = int(domkey.om)
                if domkey.pmt > self.npmts  :
                    self.npmts  = int(domkey.pmt)

            DOM_space = domsUsed[OMKey(1,1,0)].position.z - domsUsed[OMKey(1,2,0)].position.z

            String_space = 99999.
            string_pos = list()

            for i in range(1,self.nstrings+1) :
                string_pos.append(domsUsed[OMKey(i,1,0)].position)

            average_min_stringdist = 0.0
            for i in range(len(string_pos)-1):
                this_min_stringdist = 99999.
                domposi = string_pos[i]
                for j in range(i+1,len(string_pos)) :
                    domposj = string_pos[j]
                    dist = sqrt((domposi.x-domposj.x)**2.0+(domposi.y-domposj.y)**2.0)
                    if dist < this_min_stringdist :
                        this_min_stringdist = dist
                average_min_stringdist += this_min_stringdist
            average_min_stringdist /= (self.nstrings-1)


            maxDOMDistance = 0.0
            dom_pos = list()
            for domkey in domsUsed.keys() :
                dom_pos.append(domsUsed[domkey].position)

            for i in range(len(dom_pos)-1) :
                domposi = dom_pos[i]
                for j in range(i+1,len(dom_pos)):
                    domposj = dom_pos[j]
                    dist = sqrt((domposi.x-domposj.x)**2.0+(domposi.y-domposj.y)**2.0+(domposi.z-domposj.z)**2.0)
                    if dist > maxDOMDistance :
                        maxDOMDistance = dist

            self.FullDetectorCoincidenceWindow *= maxDOMDistance/0.3 + DOM_space*1.3/0.3
            self.StringCoincidenceWindow *=  average_min_stringdist/0.3 + DOM_space*1.3/0.3


        self.StringTriggerGroups = list()

        #Figure out trigger groups.
        if self.ForceAdjacency :

            nstart = int((self.StringNRows-1)/2)
            for j in range(self.nstrings) :
                for i in range(nstart,self.nDOMs-nstart) :
                    self.StringTriggerGroups.append([])
                    for l in range(len(string_pos)) :
                        if sqrt((string_pos[j].x-string_pos[l].x)**2.0+(string_pos[j].y-string_pos[l].y)**2.0) < average_min_stringdist*1.5 :
                            for k in range(i-nstart,i+nstart+1) :
                                self.StringTriggerGroups[-1].append(OMKey(l+1,k+1,0))

        else :
            nstart = int((self.StringNRows-1)/2)
            for i in range(nstart,self.nDOMs-nstart) :
                self.StringTriggerGroups.append([])
                for j in range(1,self.nstrings+1) :
                    for k in range(i-nstart,i+nstart+1) :
                        self.StringTriggerGroups[-1].append(OMKey(j,k,0))  

        self.eventcount = 0
    def DetectorStatus(self,frame) :

        frame["FullDetectorCoincidenceWindow"+self.output] = dataclasses.I3Double(self.FullDetectorCoincidenceWindow)
        frame["StringCoincidenceWindow"+self.output] = dataclasses.I3Double(self.StringCoincidenceWindow)
        frame["FullDetectorCoincidenceN"+self.output] = dataclasses.I3Double(self.FullDetectorCoincidenceN) 
        frame["StringCoincidenceN"+self.output] = dataclasses.I3Double(self.StringCoincidenceN)
        frame["TriggerForceAdjacency"+self.output] = dataclasses.I3Double(self.ForceAdjacency)
    
        self.PushFrame(frame)

    def GetDOMTriggers(self,DOMCoincidence_time,DOMCoincidence_ncoin,DOMCoincidence_pmts,ncoincidence) :

        DOMTriggers = {}

        for key in DOMCoincidence_time.keys() :
            fullpmtlist = DOMCoincidence_pmts[key]
            start=0
            for i in range(len(DOMCoincidence_time[key])) :
                time = DOMCoincidence_time[key][i]
                coinc = DOMCoincidence_ncoin[key][i] 
                if coinc >= ncoincidence :
                    if key not in DOMTriggers.keys() :
                        DOMTriggers[key] = list()
                    DOMTriggers[key].append(time)
                    #OMKEY_pmt = OMkey(key.string,key.om,i)
                    #DOMTriggers[key][time] = DOMCoincidence_pmts[OMKEY_pmt]

        return DOMTriggers
        
    
    def DAQ(self,frame) :

        DOMCoincidence_time = frame["DOMTrigger_time"+self.input]
        DOMCoincidence_ncoin = frame["DOMTrigger_ncoin"+self.input]
        DOMCoincidence_pmts = frame["DOMTrigger_pmts"+self.input]

        DOMTriggers = self.GetDOMTriggers(DOMCoincidence_time,DOMCoincidence_ncoin,DOMCoincidence_pmts,self.DOMPMTCoinc)

        #print("DOMTriggers")
        #print(DOMTriggers)

        FullDetectDOMTriggers = list()
        StringTriggers = {}
    
        for key in DOMTriggers.keys() :
            for i in range(len(self.StringTriggerGroups)) :
                if key in self.StringTriggerGroups[i]:
                    if i not in StringTriggers.keys() :
                        StringTriggers[i] = list()
                    for time in DOMTriggers[key] :
                        StringTriggers[i].append((key,time))
            FullDetectDOMTriggers.append((key,time))

    
        StringTrigOpp = {}
        if self.DoStringTrigger :
            for i in StringTriggers.keys() :
                StringTrigOpp[i] = list()
                for j in range(len(StringTriggers[i])) :
                    StringTrigOpp[i].append([StringTriggers[i][j][1],[StringTriggers[i][j][0]],[StringTriggers[i][j][1]]])
                for j in range(len(StringTriggers[i])) :
                    for k in range(len(StringTrigOpp[i])) :
                        if StringTriggers[i][j][1] - StringTrigOpp[i][k][0] < self.StringCoincidenceWindow and StringTriggers[i][j][1] - StringTrigOpp[i][k][0] >= 0.0 and StringTriggers[i][j][0] not in StringTrigOpp[i][k][1]:
                            StringTrigOpp[i][k][1].append(StringTriggers[i][j][0])
                            StringTrigOpp[i][k][2].append(StringTriggers[i][j][1])

        #print("FullDetectDOMTriggers")
        #print(FullDetectDOMTriggers)
        DetectTrigOpp = list()
        for j in range(len(FullDetectDOMTriggers)) :
            DetectTrigOpp.append([FullDetectDOMTriggers[j][1],[FullDetectDOMTriggers[j][0]],[FullDetectDOMTriggers[j][1]]])
        for j in range(len(FullDetectDOMTriggers)) :
            for k in range(len(DetectTrigOpp)) :
                if ((FullDetectDOMTriggers[j][1] - DetectTrigOpp[k][0]) < self.FullDetectorCoincidenceWindow) and ((FullDetectDOMTriggers[j][1] - DetectTrigOpp[k][0]) >= 0.0) and (FullDetectDOMTriggers[j][0] not in DetectTrigOpp[k][1]):
                        DetectTrigOpp[k][1].append(FullDetectDOMTriggers[j][0])
                        DetectTrigOpp[k][2].append(FullDetectDOMTriggers[j][1])

        #print("DetectTrigOpp")
        #print(DetectTrigOpp)

        stringTriggerTime = dataclasses.I3VectorDouble()
        detectorTriggerTime = dataclasses.I3VectorDouble()
    
        for i in range(len(DetectTrigOpp)):
            if len(DetectTrigOpp[i][1]) >= self.FullDetectorCoincidenceN :
                triggered = False
                for k in range(len(detectorTriggerTime)):
                    if abs(detectorTriggerTime[k] - max(DetectTrigOpp[i][2]))<self.EventLength :
                        detectorTriggerTime[k] = min(detectorTriggerTime[k],max(DetectTrigOpp[i][2]))
                        triggered = True
                if not triggered :
                    detectorTriggerTime.append(max(DetectTrigOpp[i][2]))

        if self.DoStringTrigger :
            for j in StringTrigOpp.keys() :
                for i in range(len(StringTrigOpp[j])):
                    if len(StringTrigOpp[j][i][1]) >= self.StringCoincidenceN :
                        triggered = False
                        for k in range(len(stringTriggerTime)):
                            if abs(stringTriggerTime[k] - max(StringTrigOpp[j][i][2]))<self.EventLength :
                                stringTriggerTime[k] = min(stringTriggerTime[k],max(StringTrigOpp[j][i][2]))
                                triggered = True
                        if not triggered :
                            stringTriggerTime.append(max(StringTrigOpp[j][i][2]))

        if self.CutOnTrigger and len(stringTriggerTime) < 1 and len(detectorTriggerTime) < 1 :
            if self.CutOnTrigger :
                return
            else :
                frame["DetectorTriggers"+self.output] = detectorTriggerTime
                frame["StringTriggers"+self.output] = stringTriggerTime
                outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
                frame[self.PulseSeriesOut] = outputpulsemap
                self.PushFrame(frame)

        #print("detectorTriggerTime")
        #print(detectorTriggerTime)
        #print("StringTriggerTime")
        #print(StringTriggerTime)

        frame["DetectorTriggers"+self.output] = detectorTriggerTime
        frame["StringTriggers"+self.output] = stringTriggerTime

        pulseseriesmap = frame[self.PulseSeriesIn]
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()

        mintrigtime = min(detectorTriggerTime,stringTriggerTime)

        for dom in pulseseriesmap :
            pulseseries = dataclasses.I3RecoPulseSeries()
            for pulse in pulseseriesmap[dom] :
                if pulse.time > mintrigtime-self.TriggerTime and pulse.time < mintrigtime+self.EventLength-self.TriggerTime :
                    resetpulse = dataclasses.I3RecoPulse()
                    resetpulse.charge = pulse.charge
                    resetpulse.time = pulse.time-mintrigtime+self.TriggerTime
                    pulseseries.append(resetpulse)
            if len(pulseseries) > 0 :
                outputseriesmap[dom] = pulseseries

        frame[self.PulseSeriesOut] = outputpulsemap
                
        self.PushFrame(frame)
