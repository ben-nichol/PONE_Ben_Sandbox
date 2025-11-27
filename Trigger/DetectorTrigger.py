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
        self.AddParameter("output", "Append the outputs", "")
        self.AddParameter("input", "Name of the Physics I3MCTree name", "")
        self.AddParameter("CutOnTrigger", "Cut events that do not trigger.", False)
        self.AddParameter("FullDetectorCoincidenceN", "", 3)
        self.AddParameter("FullDetectorCoincidenceWindow", "", 1.1)
        self.AddParameter("StringCoincidenceN", "", 2)
        self.AddParameter("StringCoincidenceWindow", "", 1.1)
        self.AddParameter("StringNRows", "", 3)
        self.AddParameter("StringDist", "", 1.5)
        self.AddParameter(
            "ScaleBySpacing",
            "Turn the Windows inputs to be relative to detector size.",
            True,
        )
        self.AddParameter("ForceAdjacency", "Require adjacency ", True)
        self.AddParameter("OMPMTCoinc", "Number of PMTs needed for OM Trigger", 2)
        self.AddParameter("EventLength", "Length of Event", 10000)
        self.AddParameter("TriggerTime", "Time of trigger in event.", 2000)
        self.AddParameter("PulseSeriesIn", "Pulse series in", "")
        self.AddParameter("PulseSeriesOut", "Pulse series out", "")
        self.AddParameter("SingleOMTriggerCoince", " ", 3)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.output = self.GetParameter("output")
        self.input = self.GetParameter("input")
        self.CutOnTrigger = self.GetParameter("CutOnTrigger")
        self.FullDetectorCoincidenceN = self.GetParameter("FullDetectorCoincidenceN")
        self.StringCoincidenceN = self.GetParameter("StringCoincidenceN")
        self.ForceAdjacency = self.GetParameter("ForceAdjacency")
        self.StringNRows = self.GetParameter("StringNRows")
        self.StringDist = self.GetParameter("StringDist")
        self.FullDetectorCoincidenceWindow_unscaled = self.GetParameter(
            "FullDetectorCoincidenceWindow"
        )
        self.StringCoincidenceWindow_unscaled = self.GetParameter(
            "StringCoincidenceWindow"
        )
        self.OMPMTCoinc = self.GetParameter("OMPMTCoinc")
        self.EventLength = self.GetParameter("EventLength")
        self.TriggerTime = self.GetParameter("TriggerTime")
        self.PulseSeriesIn = self.GetParameter("PulseSeriesIn")
        self.PulseSeriesOut = self.GetParameter("PulseSeriesOut")
        self.ScaleBySpacing = self.GetParameter("ScaleBySpacing")
        self.DoStringTrigger = True
        self.DoCoincTriggers = True
        self.SingleOMPMTTriggerCoince = self.GetParameter("SingleOMTriggerCoince")
        if self.SingleOMPMTTriggerCoince <= self.OMPMTCoinc:
            self.DoCoincTriggers = False

        self.nstrings = int(0)
        self.nOMs = int(0)
        self.npmts = int(0)

        self.has_seen_geometry = False

        self.eventcount = 0

    def Geometry(self, frame):
        self.has_seen_geometry = True

        if self.ScaleBySpacing:
            # Figure out largest distance between two OMs, average distance between closest String, Average min distance between OMs.
            geo = frame["I3Geometry"].omgeo

            stringlist = set()
            omlist = set()
            pmtlist = set()

            for omkey in geo.keys():
                stringlist.add(omkey.string)
                omlist.add(omkey.om)
                pmtlist.add(omkey.pmt)

            stringlist=sorted(stringlist)
            self.nstrings = len(stringlist)
            self.nOMs = len(omlist)
            self.npmts = len(pmtlist) #assumes all POM and no PCAL

            keys = list(geo.keys())
            OM_space = abs(geo[keys[0]].position.z- geo[keys[self.npmts]].position.z) #npmts should index to the first pmt of the next OM
            string_pos = list()

            for string in stringlist:
                string_pos.append(geo[OMKey(string, 1, 1)].position)

            average_min_stringdist = 0.0
            for i in range(len(string_pos) - 1):
                this_min_stringdist = 99999.0
                omposi = string_pos[i]
                for j in range(i + 1, len(string_pos)):
                    omposj = string_pos[j]
                    dist = sqrt((omposi.x - omposj.x) ** 2.0 + (omposi.y - omposj.y) ** 2.0)
                    if dist < this_min_stringdist:
                        this_min_stringdist = dist
                average_min_stringdist += this_min_stringdist
            average_min_stringdist /= self.nstrings - 1

            maxOMDistance = 0.0
            om_pos = list()
            for omkey in geo.keys():
                if (omkey.om==1 or omkey.om==self.nOMs): #assumes 1 to N of POMs 
                    if(omkey.pmt==1): #only check 1 PMT
                        om_pos.append(geo[omkey].position)

            for i in range(len(om_pos) - 1):
                omposi = om_pos[i]
                for j in range(i + 1, len(om_pos)):
                    omposj = om_pos[j]
                    dist = sqrt(
                          (omposi.x - omposj.x) ** 2.0
                        + (omposi.y - omposj.y) ** 2.0
                        + (omposi.z - omposj.z) ** 2.0
                    )
                    if dist > maxOMDistance:
                        maxOMDistance = dist

            self.FullDetectorCoincidenceWindow = (
                self.FullDetectorCoincidenceWindow_unscaled * maxOMDistance / 0.3
                + OM_space * 1.3 / 0.3
            )
            self.StringCoincidenceWindow = (
                self.StringCoincidenceWindow_unscaled * average_min_stringdist / 0.3
                + OM_space * 1.3 / 0.3
            )
        else:
            self.FullDetectorCoincidenceWindow = (
                self.FullDetectorCoincidenceWindow_unscaled
            )
            self.StringCoincidenceWindow = self.StringCoincidenceWindow_unscaled

        if self.StringCoincidenceN >= self.FullDetectorCoincidenceN:
            if self.StringCoincidenceWindow <= self.FullDetectorCoincidenceWindow:
                self.DoStringTrigger = False

        self.StringTriggerGroups = list()

        # Figure out trigger groups.
        if self.ForceAdjacency:
            nstart = int((self.StringNRows - 1) / 2)
            for j in range(self.nstrings):
                for i in range(nstart, self.nOMs - nstart):
                    self.StringTriggerGroups.append([])
                    for l in range(len(string_pos)):
                        if (
                            sqrt(
                                (string_pos[j].x - string_pos[l].x) ** 2.0
                                + (string_pos[j].y - string_pos[l].y) ** 2.0
                            )
                            < average_min_stringdist * 1.5
                        ):
                            for k in range(i - nstart, i + nstart + 1):
                                self.StringTriggerGroups[-1].append(OMKey(stringlist[l], k + 1, 1))
        else:
            nstart = int((self.StringNRows - 1) / 2)
            for i in range(nstart, self.nOMs - nstart):
                self.StringTriggerGroups.append([])
                for j in range(self.nstrings):
                    for k in range(i - nstart, i + nstart + 1):
                        self.StringTriggerGroups[-1].append(OMKey(stringlist[j], k, 1))

        self.OMTriggerGroups = {}
        for i in range(len(self.StringTriggerGroups)):
            for om in self.StringTriggerGroups[i]:
                if om not in self.OMTriggerGroups.keys():
                    self.OMTriggerGroups[om] = list()
                self.OMTriggerGroups[om].append(i)

        self.PushFrame(frame)

    def DetectorStatus(self, frame):
        if not self.has_seen_geometry:
            raise RuntimeError(
                "This module needs a Geometry frame in your input stream"
            )

        frame["FullDetectorCoincidenceWindow" + self.output] = dataclasses.I3Double(
            self.FullDetectorCoincidenceWindow
        )
        frame["StringCoincidenceWindow" + self.output] = dataclasses.I3Double(
            self.StringCoincidenceWindow
        )
        frame["FullDetectorCoincidenceN" + self.output] = dataclasses.I3Double(
            self.FullDetectorCoincidenceN
        )
        frame["StringCoincidenceN" + self.output] = dataclasses.I3Double(
            self.StringCoincidenceN
        )
        frame["TriggerForceAdjacency" + self.output] = dataclasses.I3Double(
            self.ForceAdjacency
        )

        self.PushFrame(frame)

    def GetOMTriggers(
        self,
        OMCoincidence_time,
        OMCoincidence_ncoin,
        OMCoincidence_pmts,
        ncoincidence,
    ):
        OMTriggers = {}

        for key in OMCoincidence_time.keys():
            fullpmtlist = OMCoincidence_pmts[key]
            start = 0
            for i in range(len(OMCoincidence_time[key])):
                time = OMCoincidence_time[key][i]
                coinc = OMCoincidence_ncoin[key][i]
                if coinc >= ncoincidence:
                    if key not in OMTriggers.keys():
                        OMTriggers[key] = list()
                    OMTriggers[key].append(time)

        return OMTriggers

    def DAQ(self, frame):
        # print("Detecctor Trigger Start")

        OMCoincidence_time = frame["DOMTrigger_Time" + self.input]
        OMCoincidence_ncoin = frame["DOMTrigger_NCoin" + self.input]
        OMCoincidence_pmts = frame["DOMTrigger_PMTs" + self.input]

        stringTriggerTime = dataclasses.I3VectorDouble()
        detectorTriggerTime = dataclasses.I3VectorDouble()
        singleOMTriggerTime = dataclasses.I3VectorDouble()

        SingleOMTriggers = self.GetOMTriggers(
            OMCoincidence_time, OMCoincidence_ncoin, OMCoincidence_pmts, 3
        )

        if self.DoCoincTriggers:
            OMTriggers = self.GetOMTriggers(
                OMCoincidence_time,
                OMCoincidence_ncoin,
                OMCoincidence_pmts,
                self.OMPMTCoinc,
            )

            FullDetectOMTriggers = list()
            StringTriggers = {}

            # for key in OMTriggers.keys() :
            #    for i in range(len(self.StringTriggerGroups)) :
            #        if key in self.StringTriggerGroups[i]:
            #            if i not in StringTriggers.keys() :
            #                StringTriggers[i] = list()
            #            for time in OMTriggers[key] :
            #                StringTriggers[i].append((key,time))
            #    for time in OMTriggers[key] :
            #        FullDetectOMTriggers.append((key,time))

            for key in OMTriggers.keys():
                for i in self.OMTriggerGroups[key]:
                    if i not in StringTriggers.keys():
                        StringTriggers[i] = list()
                    for time in OMTriggers[key]:
                        StringTriggers[i].append((key, time))
                for time in OMTriggers[key]:
                    FullDetectOMTriggers.append((key, time))

            StringTrigOpp = {}
            if self.DoStringTrigger:
                # print(StringTriggers.keys())
                for i in StringTriggers.keys():
                    StringTrigOpp[i] = list()
                    for j in range(len(StringTriggers[i])):
                        StringTrigOpp[i].append(
                            [
                                StringTriggers[i][j][1],
                                [StringTriggers[i][j][0]],
                                [StringTriggers[i][j][1]],
                            ]
                        )
                    for k in range(len(StringTrigOpp[i])):
                        for j in range(len(StringTriggers[i])):
                            if (
                                abs(StringTriggers[i][j][1] - StringTrigOpp[i][k][0])
                                < self.StringCoincidenceWindow
                                and StringTriggers[i][j][0]
                                not in StringTrigOpp[i][k][1]
                            ):
                                StringTrigOpp[i][k][1].append(StringTriggers[i][j][0])
                                StringTrigOpp[i][k][2].append(StringTriggers[i][j][1])

            # print("FullDetectOMTriggers")
            # print(FullDetectOMTriggers)
            DetectTrigOpp = list()
            for j in range(len(FullDetectOMTriggers)):
                DetectTrigOpp.append(
                    [
                        FullDetectOMTriggers[j][1],
                        [FullDetectOMTriggers[j][0]],
                        [FullDetectOMTriggers[j][1]],
                    ]
                )
            for k in range(len(DetectTrigOpp)):
                for j in range(len(FullDetectOMTriggers)):
                    if (
                        (
                            (FullDetectOMTriggers[j][1] - DetectTrigOpp[k][0])
                            < self.FullDetectorCoincidenceWindow
                        )
                        and ((FullDetectOMTriggers[j][1] - DetectTrigOpp[k][0]) >= 0.0)
                        and (FullDetectOMTriggers[j][0] not in DetectTrigOpp[k][1])
                    ):
                        DetectTrigOpp[k][1].append(FullDetectOMTriggers[j][0])
                        DetectTrigOpp[k][2].append(FullDetectOMTriggers[j][1])

            # print("DetectTrigOpp")
            # print(DetectTrigOpp)

            for i in range(len(DetectTrigOpp)):
                if len(DetectTrigOpp[i][1]) >= self.FullDetectorCoincidenceN:
                    triggered = False
                    for k in range(len(detectorTriggerTime)):
                        if (
                            abs(detectorTriggerTime[k] - max(DetectTrigOpp[i][2]))
                            < self.EventLength
                        ):
                            detectorTriggerTime[k] = min(
                                detectorTriggerTime[k], max(DetectTrigOpp[i][2])
                            )
                            triggered = True
                    if not triggered:
                        detectorTriggerTime.append(max(DetectTrigOpp[i][2]))

            if self.DoStringTrigger:
                for j in StringTrigOpp.keys():
                    for i in range(len(StringTrigOpp[j])):
                        if len(StringTrigOpp[j][i][1]) >= self.StringCoincidenceN:
                            triggered = False
                            for k in range(len(stringTriggerTime)):
                                if (
                                    abs(
                                        stringTriggerTime[k]
                                        - max(StringTrigOpp[j][i][2])
                                    )
                                    < self.EventLength
                                ):
                                    stringTriggerTime[k] = min(
                                        stringTriggerTime[k],
                                        max(StringTrigOpp[j][i][2]),
                                    )
                                    triggered = True
                            if not triggered:
                                stringTriggerTime.append(max(StringTrigOpp[j][i][2]))

        for om in SingleOMTriggers.keys():
            for time in SingleOMTriggers[om]:
                singleOMTriggerTime.append(time)

        if (
            self.CutOnTrigger
            and len(stringTriggerTime) < 1
            and len(detectorTriggerTime) < 1
            and len(singleOMTriggerTime) < 1
        ):
            if self.CutOnTrigger:
                return
            else:
                frame["DetectorTriggers" + self.output] = detectorTriggerTime
                frame["StringTriggers" + self.output] = stringTriggerTime
                frame["singleOMTrigger" + self.output] = singleOMTriggerTime
                outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
                frame[self.dPulseSeriesOut] = outputpulsemap
                self.PushFrame(frame)

        frame["DetectorTriggers" + self.output] = detectorTriggerTime
        frame["StringTriggers" + self.output] = stringTriggerTime
        frame["singleOMTrigger" + self.output] = singleOMTriggerTime

        pulseseriesmap = frame[self.PulseSeriesIn]
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()

        mintrigtime = 99999999.0
        for _time in detectorTriggerTime:
            if _time < mintrigtime:
                mintrigtime = _time

        for _time in stringTriggerTime:
            if _time < mintrigtime:
                mintrigtime = _time

        for _time in singleOMTriggerTime:
            if _time < mintrigtime:
                mintrigtime = _time

        for om in pulseseriesmap.keys():
            pulseseries = dataclasses.I3RecoPulseSeries()
            for pulse in pulseseriesmap[om]:
                if (
                    pulse.time > mintrigtime - self.TriggerTime
                    and pulse.time < mintrigtime + self.EventLength - self.TriggerTime
                ):
                    resetpulse = dataclasses.I3RecoPulse()
                    resetpulse.charge = pulse.charge
                    resetpulse.time = pulse.time - mintrigtime + self.TriggerTime
                    pulseseries.append(resetpulse)
            if len(pulseseries) > 0:
                outputpulsemap[om] = pulseseries

        frame[self.PulseSeriesOut] = outputpulsemap
        frame["TriggerTime" + self.output] = dataclasses.I3Double(mintrigtime)

        # print("Detector Trigger Done")

        self.PushFrame(frame)
        Pframe = icetray.I3Frame("P")
        self.PushFrame(Pframe)
