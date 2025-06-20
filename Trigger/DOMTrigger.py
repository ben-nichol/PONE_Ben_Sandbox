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
        self.AddParameter("inputmap", "Name of the Physics I3MCTree name", "I3RecoPulseSeriesMap")
        self.AddParameter("trigger_map", "Name of om wide pulse series map to trigger off of", "triggerpulsemap")
        self.AddParameter("PEthreshold", " Pulse charge threshold", 0.25)
        self.AddParameter("CutNotTriggered", "Cut events that do not trigger.", False)
        self.AddParameter("SingleDOMCoincidenceN", "", 3)
        self.AddParameter("SingleDOMCoincidenceWindow", "", 10)
        self.AddParameter("SingleStringNRows", "", 3)
        self.AddParameter("ForceAdjacency", "Require adjacency ", True)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.output                     = self.GetParameter("output")
        self.inputmap                   = self.GetParameter("inputmap")
        self.trigger_map                = self.GetParameter("trigger_map")
        self.PEthreshold                = self.GetParameter("PEthreshold")
        self.CutNotTriggered            = self.GetParameter('CutNotTriggered')
        self.SingleDOMCoincidenceN      = self.GetParameter("SingleDOMCoincidenceN")
        self.SingleDOMCoincidenceWindow = self.GetParameter("SingleDOMCoincidenceWindow")

        self.FRAME_NUM = 1 ########################################################################################################

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
        # print(self.FRAME_NUM) ########################################################################################################

        PulseSeriesMap = frame[self.trigger_map]
        DOMCoincidence_dict = {}
        DOMPulseCount = {}
        DOMPMTCount = {}

        DOMCoincidence_time      = dataclasses.I3MapKeyVectorDouble() # times stored in this vector are actually integers because they are used as dict keys
        DOMCoincidence_hit_times = dataclasses.I3MapKeyVectorDouble()
        DOMCoincidence_ncoin     = dataclasses.I3MapKeyVectorInt()
        DOMCoincidence_pmts      = dataclasses.I3MapKeyVectorInt()

        for omkey in PulseSeriesMap.keys():
            DOMCoincidence_dict[omkey] = {}
            DOMPMTCount[omkey]         = set()
            pulses                     = PulseSeriesMap[omkey]
            lookback                   = 0

            for i, pulse in enumerate(pulses):
                pmt         = int(pulse.width)
                key_time    = int(pulse.time)
                actual_time = pulse.time
                if key_time not in DOMCoincidence_dict[omkey].keys():
                    DOMCoincidence_dict[omkey][key_time] = {}
                    DOMCoincidence_dict[omkey][key_time]['pmts'] = {pmt}
                    DOMCoincidence_dict[omkey][key_time]['pmt_time_info'] = {(pmt, actual_time)}
                    DOMPMTCount[omkey].add(pmt)
                else:

                    DOMPMTCount[omkey].add(pmt)
                    # [pmts] is a set so we can add no problem and not get duplicates but actual times
                    # will probably always be different so it won't add a pmt but will add a time if
                    # we don't check here
                    if pmt not in DOMCoincidence_dict[omkey][key_time]['pmts']:
                        DOMCoincidence_dict[omkey][key_time]['pmts'].add(pmt)
                        DOMCoincidence_dict[omkey][key_time]['pmt_time_info'].add((pmt, actual_time))

                for j in range(lookback, i):
                    backpulse     = pulses[j]
                    backpulse_pmt = int(backpulse.width)
                    if backpulse.time < actual_time - self.SingleDOMCoincidenceWindow:
                        lookback = j
                    else:
                        if backpulse_pmt not in DOMCoincidence_dict[omkey][key_time]['pmts']:
                            DOMCoincidence_dict[omkey][key_time]['pmts'].add(backpulse_pmt)
                            DOMCoincidence_dict[omkey][key_time]['pmt_time_info'].add((backpulse_pmt, backpulse.time))

            if len(DOMPMTCount[omkey]) >= self.SingleDOMCoincidenceN:
                times     = list()
                hit_times = list()
                coinc     = list()
                pmts      = list()
                for time in DOMCoincidence_dict[omkey].keys():
                    if len(DOMCoincidence_dict[omkey][time]['pmts']) < self.SingleDOMCoincidenceN:
                        continue
                    times.append(time)
                    coinc.append(len(DOMCoincidence_dict[omkey][time]['pmts']))
                    for pmt, hit_time in DOMCoincidence_dict[omkey][time]['pmt_time_info']:
                        pmts.append(int(pmt))
                        hit_times.append(hit_time)
                if len(times) > 0:
                    DOMCoincidence_time[omkey]      = times
                    DOMCoincidence_hit_times[omkey] = hit_times
                    DOMCoincidence_ncoin[omkey]     = coinc
                    DOMCoincidence_pmts[omkey]      = pmts

                # print(times)
        if self.CutNotTriggered:
            if len(DOMCoincidence_ncoin) == 0:
                self.FRAME_NUM += 1 ########################################################################################################
                return


        frame["DOMTrigger_Time" + self.output]      = DOMCoincidence_time
        frame["DOMTrigger_NCoin" + self.output]     = DOMCoincidence_ncoin
        frame["DOMTrigger_PMTs" + self.output]      = DOMCoincidence_pmts
        frame["DOMTrigger_HitTimes" + self.output] = DOMCoincidence_hit_times
        # print("DOM Trigger End")

        self.PushFrame(frame)

        self.FRAME_NUM += 1 ########################################################################################################
