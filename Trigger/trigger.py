from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from math import sqrt
from copy import deepcopy

class Trigger(icetray.I3ConditionalModule):
  """
  Simple Implementation of the PMT response.
  """

  def __init__(self, context):
    icetray.I3ConditionalModule.__init__(self, context)
    self.AddParameter("GCDFile","GCD to be simulated",'')
    self.AddParameter("output","Append the outputs",'')
    self.AddParameter("inputmap","Name of the Physics I3MCTree name","I3RecoPulseSeriesMap")
    self.AddParameter("PEthreshold"," Pulse charge threshold",0.25)
    self.AddParameter("CutOnTrigger","Cut events that do not trigger.",False)
    self.AddParameter("FullDetectorCoincidenceN","",3)
    self.AddParameter("FullDetectorCoincidenceWindow","",1.1)
    self.AddParameter("StringCoincidenceN","",2)
    self.AddParameter("StringCoincidenceWindow","",1.1)
    self.AddParameter("InterStringCoincidenceN","",2)
    self.AddParameter("InterStringCoincidenceWindow","",1.1)
    self.AddParameter("InterStringNRows","",3)
    self.AddParameter("InterStringDist","",1.5)
    self.AddParameter("SingleDOMCoincidenceN","",2)
    self.AddParameter("SingleDOMCoincidenceWindow","",10)
    self.AddParameter("SingleStringNRows","",3)
    self.AddParameter("OMPlaneCoincidenceN","",2)
    self.AddParameter("OMPlaneCoincidenceWindow","",1.1)
    self.AddParameter("ScaleBySpacing","Turn the Windows inputs to be relative to detector size.",True)
    self.AddParameter("ForceAdjacency","Require adjacency ",True)
    self.AddOutBox("OutBox")

  def Configure(self):

    self.gcdFile = self.GetParameter("GCDFile")
    self.output = self.GetParameter("output")
    self.inputmap = self.GetParameter("inputmap")
    self.PEthreshold = self.GetParameter("PEthreshold")
    self.CutOnTrigger = self.GetParameter("CutOnTrigger")
    self.FullDetectorCoincidenceN = self.GetParameter("FullDetectorCoincidenceN")
    self.StringCoincidenceN = self.GetParameter("StringCoincidenceN")
    self.InterStringCoincidenceN = self.GetParameter("InterStringCoincidenceN")
    self.SingleDOMCoincidenceN = self.GetParameter("SingleDOMCoincidenceN")
    self.SingleDOMCoincidenceWindow = self.GetParameter("SingleDOMCoincidenceWindow")
    self.OMPlaneCoincidenceN = self.GetParameter("OMPlaneCoincidenceN")
    self.OMPlaneCoincidenceWindow = self.GetParameter("OMPlaneCoincidenceWindow")
    self.ForceAdjacency = self.GetParameter("ForceAdjacency")
    self.SingleStringNRows = self.GetParameter("SingleStringNRows")
    self.InterStringNRows = self.GetParameter("InterStringNRows")
    self.InterStringDist = self.GetParameter("InterStringDist")
    self.FullDetectorCoincidenceWindow = self.GetParameter("FullDetectorCoincidenceWindow")
    self.InterStringCoincidenceWindow = self.GetParameter("InterStringCoincidenceWindow")
    self.StringCoincidenceWindow = self.GetParameter("StringCoincidenceWindow")

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

      self.nstrings += 1
      self.nDOMs += 1
      self.npmts += 1
      
      DOM_space = domsUsed[OMKey(0,0,0)].position.z - domsUsed[OMKey(0,1,0)].position.z

      String_space = 99999.
      string_pos = list()

      for i in range(self.nstrings) :
        string_pos.append(domsUsed[OMKey(i,0,0)].position)

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
      self.InterStringCoincidenceWindow *=  average_min_stringdist/0.3 + DOM_space*1.3/0.3
      self.StringCoincidenceWindow *= 2.0*DOM_space*1.3/0.3


    self.SingleStringTriggerGroups = list()
    self.InterStringTriggerGroups = list()

    #Figure out trigger groups.
    if self.ForceAdjacency :

      nstart = int((self.SingleStringNRows-1)/2)
      for j in range(self.nstrings) :
        for i in range(nstart,self.nDOMs-nstart) :
          self.SingleStringTriggerGroups.append([])
          for l in range(i-nstart,i+nstart+1) :
            self.SingleStringTriggerGroups[-1].append(OMKey(j,l,0))

      nstart = int((self.InterStringNRows-1)/2)
      for j in range(self.nstrings) :
        for i in range(nstart,self.nDOMs-nstart) :
          self.InterStringTriggerGroups.append([])
          for l in range(len(string_pos)) :
            if sqrt((string_pos[j].x-string_pos[l].x)**2.0+(string_pos[j].y-string_pos[l].y)**2.0) < average_min_stringdist*1.5 :
              for k in range(i-nstart,i+nstart+1) :
                self.InterStringTriggerGroups[-1].append(OMKey(l,k,0))

    else :
      nstart = int((self.SingleStringNRows-1)/2)
      self.SingleStringTriggerGroups.append([])
      for j in range(self.nstrings) :
        for i in range(self.nDOMs) :
          self.SingleStringTriggerGroups[-1].append(OMKey(j,i,0))

      nstart = int((self.InterStringNRows-1)/2)
      for i in range(nstart,self.nDOMs-nstart) :
        self.InterStringTriggerGroups.append([])
        for j in range(self.nstrings) :
          for k in range(i-nstart,i+nstart+1) :
            self.InterStringTriggerGroups[-1].append(OMKey(j,k,0))  


    self.Empty_FullDetectCoincidence = list()
    self.Empty_DOMCoincidence = list()
    self.Empty_StringCoincidence = list()
    self.Empty_InterStringCoincidence = list()

    for l in range(self.nstrings) :
      self.Empty_DOMCoincidence.append([])
      for k in range(self.nDOMs) :
        self.Empty_DOMCoincidence[-1].append([])
    for j in range(len(self.SingleStringTriggerGroups)) :
        self.Empty_StringCoincidence.append([])
    for j in range(len(self.InterStringTriggerGroups)) :
        self.Empty_InterStringCoincidence.append([])

    for i in range(10000) :
      self.Empty_FullDetectCoincidence.append(0.0)
      for j in range(len(self.SingleStringTriggerGroups)) :
        self.Empty_StringCoincidence[j].append(0.0)
      for j in range(len(self.InterStringTriggerGroups)) :
        self.Empty_InterStringCoincidence[j].append(0.0)
      for j in range(self.nstrings) :
        for k in range(self.nDOMs) :
          self.Empty_DOMCoincidence[j][k].append(0.0)

    self.eventcount = 0
  def DetectorStatus(self,frame) :

    frame["FullDetectorCoincidenceWindow"+self.output] = dataclasses.I3Double(self.FullDetectorCoincidenceWindow)
    frame["InterStringCoincidenceWindow"+self.output] = dataclasses.I3Double(self.InterStringCoincidenceWindow)
    frame["OMPlaneCoincidenceWindow"+self.output] = dataclasses.I3Double(self.OMPlaneCoincidenceWindow) 
    frame["FullDetectorCoincidenceN"+self.output] = dataclasses.I3Double(self.FullDetectorCoincidenceN) 
    frame["StringCoincidenceN"+self.output] = dataclasses.I3Double(self.StringCoincidenceN) 
    frame["InterStringCoincidenceN"+self.output] = dataclasses.I3Double(self.InterStringCoincidenceN)
    frame["SingleDOMCoincidenceN"+self.output] = dataclasses.I3Double(self.SingleDOMCoincidenceN)
    frame["SingleDOMCoincidenceWindow"+self.output] = dataclasses.I3Double(self.SingleDOMCoincidenceWindow)
    frame["OMPlaneCoincidenceN"+self.output] = dataclasses.I3Double(self.OMPlaneCoincidenceN)
    frame["TriggerForceAdjacency"+self.output] = dataclasses.I3Double(self.ForceAdjacency)
    
    self.PushFrame(frame)
    
  def DAQ(self,frame) :

    PulseSeriesMap = frame[self.inputmap]
    FullDetectCoincidence = deepcopy(self.Empty_FullDetectCoincidence)
    DOMCoincidence = deepcopy(self.Empty_DOMCoincidence)
    StringCoincidence = deepcopy(self.Empty_StringCoincidence)
    InterStringCoincidence = deepcopy(self.Empty_InterStringCoincidence)

    stringTrigger = False
    last_stringTriggerTime = 0
    stringTriggerTime = dataclasses.I3VectorInt()
    interStringTrigger = False
    last_interStringTriggerTime = 0
    interStringTriggerTime = dataclasses.I3VectorInt()
    detectorTrigger = False
    last_detectorTriggerTime = 0
    detectorTriggerTime = dataclasses.I3VectorInt()

    self.eventcount += 1

    npulses = 0
    for omkey in PulseSeriesMap.keys() :                                                                                                                 
        for pulse in PulseSeriesMap[omkey]: 
             npulses += 1.0    

    if npulses < self.SingleDOMCoincidenceN*min([self.StringCoincidenceN,self.InterStringCoincidenceN,self.FullDetectorCoincidenceN]) :
        if self.CutOnTrigger :
            return
        frame["DetectorTriggers"+self.output] = detectorTriggerTime
        frame["SingleStringTriggers"+self.output] = stringTriggerTime
        frame["InterStringTriggers"+self.output] = interStringTriggerTime

        self.PushFrame(frame)
        return

    for omkey in PulseSeriesMap.keys() :
      for pulse in PulseSeriesMap[omkey]:
        #if pulse.charge < 0.25 :
        #  continue

        mintimebin = max(0,min(int(pulse.time),10000))
        dommaxtimebin = max(0,min(int((pulse.time+self.SingleDOMCoincidenceWindow)),10000))

        for i in range(mintimebin,dommaxtimebin) :
          DOMCoincidence[omkey.string][omkey.om][i] += 1

    for i in range(self.nstrings) :
        for j in range(self.nDOMs) :
            intrigger=False
            omkey = OMKey(i,j,0)
            for k in range(len(DOMCoincidence[i][j])) :
                if DOMCoincidence[i][j][k] > self.SingleDOMCoincidenceN :
                    for l in range(len(self.SingleStringTriggerGroups)) :
                        if omkey in self.SingleStringTriggerGroups[l] :
                            for n in range(k,min(10000,k+int(self.StringCoincidenceWindow))):
                                StringCoincidence[l][n] += 1
                    for l in range(len(self.InterStringTriggerGroups)) :
                        if omkey in self.InterStringTriggerGroups[l] :
                            for n in range(k,min(10000,k+int(self.InterStringCoincidenceWindow))) :
                                InterStringCoincidence[l][n] += 1
                    for n in range(k,min(10000,k+int(self.FullDetectorCoincidenceWindow))) :
                        FullDetectCoincidence[n] += 1

    for i in range(10000) :
      if (i - last_stringTriggerTime) > 10000 :
        stringTrigger = False
      if not stringTrigger :
        for j in range(len(StringCoincidence)) :
          if StringCoincidence[j][i] >= self.StringCoincidenceN :
            stringTrigger = True
            last_stringTriggerTime = i
            stringTriggerTime.append(i)
            break
      if (i - last_interStringTriggerTime) > 10000 :
        interStringTrigger = False
      if not interStringTrigger :
        for j in range(len(InterStringCoincidence)) :
          if InterStringCoincidence[j][i] >= self.InterStringCoincidenceN :
            interStringTrigger = True
            last_interStringTriggerTime = i
            interStringTriggerTime.append(i)
            break
      if (i - last_detectorTriggerTime) > 10000 :
        detectorTrigger = False
      if not detectorTrigger and FullDetectCoincidence[i] >= self.FullDetectorCoincidenceN :
        last_detectorTriggerTime = i
        detectorTrigger = True
        detectorTriggerTime.append(i)

    if self.CutOnTrigger and len(stringTriggerTime) < 1 and len(detectorTriggerTime) < 1 and len(interStringTriggerTime) < 1 :
      return

    #print("detector")
    #print(detectorTriggerTime)
    #print("string")
    #print(stringTriggerTime)
    #print("interstring")
    #print(interStringTriggerTime)

    frame["DetectorTriggers"+self.output] = detectorTriggerTime
    frame["SingleStringTriggers"+self.output] = stringTriggerTime
    frame["InterStringTriggers"+self.output] = interStringTriggerTime

    self.PushFrame(frame)
