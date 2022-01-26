from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from math import sqrt
from copy import deepcopy

class DOMTrigger(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("output","Append the outputs",'')
        self.AddParameter("inputmap","Name of the Physics I3MCTree name","I3RecoPulseSeriesMap")
        self.AddParameter("PEthreshold"," Pulse charge threshold",0.25)
        self.AddParameter("CutOnTrigger","Cut events that do not trigger.",False)
        self.AddParameter("SingleDOMCoincidenceN","",2)
        self.AddParameter("SingleDOMCoincidenceWindow","",10)
        self.AddParameter("SingleStringNRows","",3)
        self.AddParameter("ForceAdjacency","Require adjacency ",True)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.output = self.GetParameter("output")
        self.inputmap = self.GetParameter("inputmap")
        self.PEthreshold = self.GetParameter("PEthreshold")
        self.SingleDOMCoincidenceN = self.GetParameter("SingleDOMCoincidenceN")
        self.SingleDOMCoincidenceWindow = self.GetParameter("SingleDOMCoincidenceWindow")

    def DetectorStatus(self,frame) :

        frame["SingleDOMCoincidenceN"+self.output] = dataclasses.I3Double(self.SingleDOMCoincidenceN)
        frame["SingleDOMCoincidenceWindow"+self.output] = dataclasses.I3Double(self.SingleDOMCoincidenceWindow)
        self.PushFrame(frame)
    
    def DAQ(self,frame) :

        PulseSeriesMap = frame[self.inputmap]
        DOMCoincidence_dict = {}
        DOMCoincidence_time = dataclasses.I3MapKeyVectorDouble()
        DOMCoincidence_ncoin = dataclasses.I3MapKeyVectorInt()
        DOMCoincidence_pmts = dataclasses.I3MapKeyVectorInt()

        #print(PulseSeriesMap)

        for omkey in PulseSeriesMap.keys() :
            OMKEY = OMKey(omkey.string,omkey.om,0)
            if OMKEY not in DOMCoincidence_dict.keys() :
                DOMCoincidence_dict[OMKEY] = {}
            for pulse in PulseSeriesMap[omkey]:
                if pulse.charge < self.PEthreshold :
                    continue
                #this is not great, what if two pulses at same time in different PMTs in same DOM???
                if float(pulse.time) not in DOMCoincidence_dict[OMKEY].keys() :
                    DOMCoincidence_dict[OMKEY][pulse.time]=[omkey.pmt]
                
                for time in DOMCoincidence_dict[OMKEY].keys() :
                    if pulse.time < float(time)+self.SingleDOMCoincidenceWindow and float(time) < pulse.time:
                        if omkey.pmt not in DOMCoincidence_dict[OMKEY][time] :
                            DOMCoincidence_dict[OMKEY][time].append(omkey.pmt)

        for OMKEY in DOMCoincidence_dict.keys() :
            times = list()
            coinc = list()
            pmts = list()
            for time in DOMCoincidence_dict[OMKEY].keys() :
                if len(DOMCoincidence_dict[OMKEY][time]) < self.SingleDOMCoincidenceN :
                    continue
                times.append(time)
                coinc.append(len(DOMCoincidence_dict[OMKEY][time]))
                for pmt in DOMCoincidence_dict[OMKEY][time] :
                    pmts.append(pmt)
            if len(times) > 0 :
                DOMCoincidence_time[OMKEY] = times
                DOMCoincidence_ncoin[OMKEY] = coinc
                DOMCoincidence_pmts[OMKEY] = pmts

        frame["DOMTrigger_time"+self.output] = DOMCoincidence_time
        frame["DOMTrigger_ncoin"+self.output] = DOMCoincidence_ncoin
        frame["DOMTrigger_pmts"+self.output] = DOMCoincidence_pmts

        self.PushFrame(frame)
