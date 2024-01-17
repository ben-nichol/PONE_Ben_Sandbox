from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars


class timeShift(icetray.I3ConditionalModule):
    """
    Shifting the timestamps of I3MCPEs to start at 7200 ns
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter(
            "MergedMCPETreeName", "Name of the Merged MCPE tree name", "MergedSeriesMap"
        )
        self.AddParameter(
            "TimeShiftedMCPE",
            "Name of the I3MCTree containing time shifted MCPEs starting at 7200 ns",
            "TimeShiftedMCPEMap",
        )
        self.AddParameter("MinTime", "Time offset for the pulses.", 7200)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.mergedSeriesName = self.GetParameter("MergedMCPETreeName")
        self.tShiftSeriesName = self.GetParameter("TimeShiftedMCPE")
        self.mintime = self.GetParameter("MinTime")

    def DAQ(self, frame):
        mcpeMap = frame[self.mergedSeriesName]
        mcpeOMKeys = mcpeMap.keys()
        timeShiftedMap = simclasses.I3CompressedPhotonSeriesMap()

        timeList = []

        for omkey in mcpeOMKeys:
            mcpeList = mcpeMap[omkey]
            timeList.extend([mcpe.time for mcpe in mcpeList])

        min_time = 0.0
        if len(timeList) != 0:
            min_time = min(timeList)
        else:
            frame[self.tShiftSeriesName] = mcpeMap
            return

        for omkey in mcpeOMKeys:
            mcpeList = mcpeMap[omkey]
            newMCPEList = simclasses.I3CompressedPhotonSeries()

            for mcpe in mcpeList:
                mcpe.time = (
                    mcpe.time - min_time
                ) + self.mintime * I3Units.ns  # [Units: ns]
                newMCPEList.append(mcpe)

            timeShiftedMap[omkey] = newMCPEList

        frame[self.tShiftSeriesName] = timeShiftedMap
        frame[self.tShiftSeriesName + str("_toffset")] = dataclasses.I3Double(
            min_time - self.mintime * I3Units.ns
        )
        self.PushFrame(frame)
