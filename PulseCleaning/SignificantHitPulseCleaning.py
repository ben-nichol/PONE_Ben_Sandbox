#!/usr/bin/env python

"""

"""

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from numpy import linalg as la
from icecube.phys_services import I3Calculator
import operator


class SignificantHitPulseCleaning(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile", "GCD to be simulated", "")
        self.AddParameter("inputseries", "GCD to be simulated", "MCPESeriesMap")
        self.AddParameter("output", "GCD to be simulated", "linefit")
        self.AddParameter("window", "Window used to decide most sig window.", 200)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.gcdFile = self.GetParameter("GCDFile")
        self.input = self.GetParameter("inputseries")
        self.output = self.GetParameter("output")
        # self.geometry = self.gcdFile.pop_frame()["I3Geometry"]
        # self.domsUsed = self.geometry.omgeo.keys()
        self.window = self.GetParameter("window")

    # A modified binary search algorithm to find all hits within a 20ns time
    # window. Takes entire list of MCPE hits, finds a 20ns window where the
    # number of hits surpasses the hit threshold, and outputs those hits. If
    # there are multiple windows, only the first occurence (in time) is taken.
    # NOTE: function assumes the list has already passed through the passDOM
    #       function. If the resulting list does not meet that criteria, the
    #       function raises a ValueError as it could cause issues for future
    #       functions relying on the significant MCPEs
    #
    # @Param:
    # mcpeList:         The list of hits the DOM registered
    # hitThresh:        The number of hits required in the window for the DOM to
    #                   be passed
    #
    # @Return:
    # An MCPESeries object with all the hits that were within the time window
    def getSignificantMCPEs(self, mcpeList):
        significantMCPEList = dataclasses.I3RecoPulseSeries()

        window_integral = {}

        npulse = 0
        for pulse in mcpeList:
            npulse += 1

        if npulse < 5:
            return mcpeList

        for pulse in mcpeList:
            time_index = pulse.time / self.window
            time_index = time_index * 5
            time_index = int(time_index)

            for i in range(time_index - 2, time_index + 3):
                if i in window_integral:
                    window_integral[i] = window_integral[i] + pulse.charge
                else:
                    window_integral[i] = pulse.charge

        max_index = int(max(window_integral.items(), key=operator.itemgetter(1))[0])

        time_index = self.window * float(max_index) / 5.0
        for pulse in mcpeList:
            if pulse.time > time_index - 2000.0 and pulse.time < time_index + 2000.0:
                significantMCPEList.append(pulse)

        return significantMCPEList

    # Converts the MCPESeriesMap to one that contains siginficant hits. Significant
    # hits are defined to be ones where the DOM passed the passDOM function. The
    # function then writes the new MCPESeriesMap to the frame before returning it
    # NOTE: function assumes the list has already passed through the passFrame
    #       function. If the resulting map does not meet that criteria, the a
    #       ValueError is raised as it could cause issues for future functions
    #       relying on the significant MCPESeriesMap
    #
    # @Param:
    # frame:            The frame to be appended then returned
    # domsUsed:         The omkeys used for the analysis. Allows to look at
    #                   smaller geometries from a larger geometry sim file
    # hitThresh:        Hit threshold for a single DOM
    # domThresh:        Number of DOMs passing the hit threshold to pass the
    #                   frame
    # maxResidual:      Maximum time residual allowed for a hit to be considered
    # GeoMap:           An I3OMGeoMap object that maps the omkeys in the
    #                   geometries to their I3OMGeo objects
    #
    # @Return:
    # The frame inputted to the function after the significant MCPESeriesMap object
    # was appended to it with the key word "MCPESerieMap_significant_hits"
    def DAQ(self, frame):
        # print("sig hit")
        mcpeMap = frame[self.input]
        significantMCPEMap = dataclasses.I3RecoPulseSeriesMap()
        # print("this is working")
        for omkey in mcpeMap.keys():
            # if omkey not in mcpeMap:
            #  continue
            if len(mcpeMap[omkey]) < 1:
                continue
            significantMCPEMap[omkey] = self.getSignificantMCPEs(mcpeMap[omkey])

        frame[self.output] = significantMCPEMap
        self.PushFrame(frame)
