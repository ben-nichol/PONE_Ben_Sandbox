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
        self.AddParameter("CutNotTriggered", "Cut events that do not trigger.", False)
        self.AddParameter("SingleDOMCoincidenceN", "", 3)
        self.AddParameter("SingleDOMCoincidenceWindow", "", 10)
        self.AddParameter("SingleStringNRows", "", 3)
        self.AddParameter("ForceAdjacency", "Require adjacency ", True)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.output = self.GetParameter("output")
        self.inputmap = self.GetParameter("inputmap")
        self.PEthreshold = self.GetParameter("PEthreshold")
        self.CutNotTriggered = self.GetParameter('CutNotTriggered')
        self.SingleDOMCoincidenceN = self.GetParameter("SingleDOMCoincidenceN")
        self.SingleDOMCoincidenceWindow = self.GetParameter(
            "SingleDOMCoincidenceWindow"
        )

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

        PulseSeriesMap = frame["triggerpulsemap"]
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
                pmt         = int(pulse.width)
                key_time    = int(pulse.time)
                actual_time = pulse.time
                if key_time not in DOMCoincidence_dict[omkey].keys():
                    DOMCoincidence_dict[omkey][key_time] = {}
                    DOMCoincidence_dict[omkey][key_time]['pmts'] = {pmt}
                    DOMCoincidence_dict[omkey][key_time]['pmt_time_info'] = {(pmt, actual_time)}
                    DOMPMTCount[omkey].add(pmt)
                else:
                    # print(self.FRAME_NUM) ########################################################################################################
                    # print('PULSE TIME ALREADY IN DICT')
                    # print(key_time)
                    # print(actual_time)
                    # print(pulse.width)
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
                        # DOMCoincidence_dict[omkey][key_time].add((int(backpulse.width), backpulse.time))
            # print(DOMPMTCount[omkey])
            if len(DOMPMTCount[omkey]) >= self.SingleDOMCoincidenceN:
                times     = list()
                hit_times = list()
                coinc     = list()
                pmts      = list()
                for time in DOMCoincidence_dict[omkey].keys():
                    if len(DOMCoincidence_dict[omkey][time]['pmts']) < self.SingleDOMCoincidenceN:
                        # print('NO COINCIDENCE') #################################################
                        continue
                    times.append(time) ################################################# put here instead of commented out below...
                    coinc.append(len(DOMCoincidence_dict[omkey][time]['pmts']))
                    for pmt, hit_time in DOMCoincidence_dict[omkey][time]['pmt_time_info']:
                        # print(pmt)
                        # print(hit_time)
                        pmts.append(int(pmt))
                        hit_times.append(hit_time)
                        #times.append(time)####################################### we need this to be the actual hit time...
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
            # total_triggers = 0
            # for omkey in DOMCoincidence_ncoin:
            #     total_triggers += len(DOMCoincidence_ncoin[omkey])
                
            # if total_triggers == 0:
            #     return


        frame["DOMTrigger_time" + self.output]      = DOMCoincidence_time
        frame["DOMTrigger_ncoin" + self.output]     = DOMCoincidence_ncoin
        frame["DOMTrigger_pmts" + self.output]      = DOMCoincidence_pmts
        frame["DOMTrigger_hit_times" + self.output] = DOMCoincidence_hit_times
        # print("DOM Trigger End")

        self.PushFrame(frame)

        self.FRAME_NUM += 1 ########################################################################################################
