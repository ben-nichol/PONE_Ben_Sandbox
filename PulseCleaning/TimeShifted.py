from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars

class timeShift(icetray.I3ConditionalModule):
    """
    Dhiting the timestamps of I3MCPEs to start at 7200 ns
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("MergedMCPETreeName",
                          "Name of the Merged MCPE tree name",
                          "MergedSeriesMap")
        self.AddParameter("TimeShiftedMCPE",
                          "Name of the I3MCTree containing time shifted MCPEs starting at 7200 ns",
                          "TimeShiftedMCPEMap")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.mergedSeriesName = self.GetParameter("MergedMCPETreeName")
        self.tShiftSeriesName = self.GetParameter("TimeShiftedMCPE")

    def DAQ(self, frame):

        mcpeMap = frame[self.mergedSeriesName]
        mcpeOMKeys = mcpeMap.keys()
        timeShiftedMap = simclasses.I3MCPESeriesMap()

        for omkey in mcpeOMKeys:
            newMCPEList = simclasses.I3MCPESeries()
            mcpeList = mcpeMap[omkey]
            timeList = np.array([mcpe.time for mcpe in mcpeList])

            if len(timeList) != 0:
                min_time = min(timeList)
                for mcpe in mcpeList:
                    mcpe.time = (mcpe.time - min_time) + 7200 #[Units: ns]
                    newMCPEList.append(mcpe)

                timeShiftedMap[omkey] = newMCPEList

        frame[self.tShiftSeriesName] = timeShiftedMap

        self.PushFrame(frame)
