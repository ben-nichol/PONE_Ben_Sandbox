from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from math import sqrt
from copy import deepcopy
import collections
from time import process_time


class DOMTrigger(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("output", "Append the outputs", "")
        self.AddParameter(
            "inputmap", "Name of the Physics I3MCTree name", "I3RecoPulseSeriesMap"
        )
        self.AddParameter("PEthreshold", " Pulse charge threshold", 0.25)
        self.AddParameter("CutOnTrigger", "Cut events that do not trigger.", False)
        self.AddParameter("SingleDOMCoincidenceN", "", 3)
        self.AddParameter("SingleDOMCoincidenceWindow", "", 10)
        self.AddParameter("SingleStringNRows", "", 3)
        self.AddParameter("ForceAdjacency", "Require adjacency ", True)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.output = self.GetParameter("output")
        self.inputmap = self.GetParameter("inputmap")
        self.PEthreshold = self.GetParameter("PEthreshold")
        self.SingleDOMCoincidenceN = self.GetParameter("SingleDOMCoincidenceN")
        self.SingleDOMCoincidenceWindow = self.GetParameter(
            "SingleDOMCoincidenceWindow"
        )

    def Geometry(self, frame):
        self.domsUsed = frame["I3Geometry"].omgeo
        self.PushFrame(frame)

    def Simulation(self, frame):
        frame["SingleDOMCoincidenceN" + self.output] = dataclasses.I3Double(
            self.SingleDOMCoincidenceN
        )
        frame["SingleDOMCoincidenceWindow" + self.output] = dataclasses.I3Double(
            self.SingleDOMCoincidenceWindow
        )
        self.PushFrame(frame)

    def DAQ(self, frame):
        # print("DOM Trigger Start")

        PulseSeriesMap = frame["triggerpulsemap"]
        DOMCoincidence_dict = {}
        DOMPulseCount = {}
        DOMPMTCount = {}

        DOMCoincidence_time = dataclasses.I3MapKeyVectorDouble()
        DOMCoincidence_ncoin = dataclasses.I3MapKeyVectorInt()
        DOMCoincidence_pmts = dataclasses.I3MapKeyVectorInt()

        for omkey in PulseSeriesMap.keys():
            DOMCoincidence_dict[omkey] = {}
            DOMPMTCount[omkey] = set()
            pulses = PulseSeriesMap[omkey]
            lookback = 0
            # if len(pulses) > 12 :
            #    averagetime = 0.0
            #    for i, pulse in enumerate(pulses) :
            #        averagetime += pulse.time
            #    averagetime /= len(pulses)
            #    minbunchtime = 9999999.0
            #    maxbunchtime = 0.0
            #    for i,pulse in enumerate(pulses) :
            #        if abs(pulse.time-averagetime)<300. :
            #            minbunchtime  = min(minbunchtime,pulse.time)
            #            maxbunchtime = max(maxbunchtime,pulse.time)
            #
            #        DOMCoincidence_time[omkey] = [minbunchtime,maxbunchtime]
            #        DOMCoincidence_ncoin[omkey] = [4,4]
            #        DOMCoincidence_pmts[omkey] = [4,4]
            #    continue
            # print("looping " + str(len(pulses)))
            for i, pulse in enumerate(pulses):
                # print(str(i) +" "+str(lookback)+" "+str(len(pulses)))
                time = int(pulse.time)
                if time not in DOMCoincidence_dict[omkey].keys():
                    DOMCoincidence_dict[omkey][time] = {int(pulse.width)}
                    DOMPMTCount[omkey].add(int(pulse.width))
                else:
                    DOMCoincidence_dict[omkey][time].add(int(pulse.width))
                    DOMPMTCount[omkey].add(int(pulse.width))

                for j in range(lookback, i):
                    backpulse = pulses[j]
                    if backpulse.time < time - self.SingleDOMCoincidenceWindow:
                        lookback = j
                    else:
                        DOMCoincidence_dict[omkey][time].add(backpulse.width)
            # print(DOMPMTCount[omkey])
            if len(DOMPMTCount[omkey]) >= self.SingleDOMCoincidenceN:
                times = list()
                coinc = list()
                pmts = list()
                for time in DOMCoincidence_dict[omkey].keys():
                    if (
                        len(DOMCoincidence_dict[omkey][time])
                        < self.SingleDOMCoincidenceN
                    ):
                        continue
                    times.append(time)
                    coinc.append(len(DOMCoincidence_dict[omkey][time]))
                    for pmt in DOMCoincidence_dict[omkey][time]:
                        pmts.append(int(pmt))
                if len(times) > 0:
                    DOMCoincidence_time[omkey] = times
                    DOMCoincidence_ncoin[omkey] = coinc
                    DOMCoincidence_pmts[omkey] = pmts

                # print(times)

        frame["DOMTrigger_time" + self.output] = DOMCoincidence_time
        frame["DOMTrigger_ncoin" + self.output] = DOMCoincidence_ncoin
        frame["DOMTrigger_pmts" + self.output] = DOMCoincidence_pmts
        # print("DOM Trigger End")

        self.PushFrame(frame)
