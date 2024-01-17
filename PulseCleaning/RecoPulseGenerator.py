from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars


class merge_recoPulse(icetray.I3ConditionalModule):
    """
    Merging MCPEs within 3 nanoseconds. Creating I3RecoPulseSeriesMap from MCPEs.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter(
            "TimeShiftedMCPE", "Name of the Physics I3MCTree name", "timeShiftedMCPEMap"
        )
        self.AddParameter(
            "RecoPulseTreeName", "Name of the noise I3MCTree name", "I3RecoPulses"
        )
        self.AddParameter(
            "MergeTime",
            "Time window within which the MCPE hits will be merged [Unit:nanoseconds]",
            3,
        )
        self.AddOutBox("OutBox")

    def Configure(self):
        self.tShiftSeriesName = self.GetParameter("TimeShiftedMCPE")
        self.recoPulseTreeName = self.GetParameter("RecoPulseTreeName")
        self.mergeTime = self.GetParameter("MergeTime")

    def DAQ(self, frame):
        timeWindow = self.mergeTime
        mcpeMap = frame[self.tShiftSeriesName]
        mcpeOMKeys = mcpeMap.keys()

        recoPulseMap = dataclasses.I3RecoPulseSeriesMap()

        for omkey in mcpeOMKeys:
            mcpeList = mcpeMap[omkey]
            timeList = [mcpe.time for mcpe in mcpeList]
            charge = np.array([mcpe.npe for mcpe in mcpeList])

            sortedCharge = [
                x for _, x in sorted(zip(timeList, charge))
            ]  # Sorting the charge with respect to the timeList
            sortedTimeList = np.sort(timeList)  # Timestamps are not sorted in timeList
            sortedmcpeList = [x for _, x in sorted(zip(timeList, mcpeList))]

            recoPulses = dataclasses.I3RecoPulseSeries()

            i = 0
            while i < len(sortedTimeList):
                mcpe = sortedmcpeList[i]
                time_end = mcpe.time + timeWindow
                times = []
                charges = []

                while mcpe.time < time_end and i < len(sortedTimeList):
                    times.append(mcpe.time)
                    charges.append(sortedCharge[i])
                    i += 1
                    if i < len(sortedTimeList):
                        mcpe = sortedmcpeList[i]

                recoPulse = dataclasses.I3RecoPulse()
                recoPulse.charge = sum(charges)
                recoPulse.width = 3 * I3Units.ns
                recoPulse.time = times[0]

                recoPulses.append(recoPulse)

            recoPulseMap[omkey] = recoPulses

        frame[self.recoPulseTreeName] = recoPulseMap

        self.PushFrame(frame)
