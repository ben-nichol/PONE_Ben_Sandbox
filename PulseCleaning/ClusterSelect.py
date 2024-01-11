# /usr/bin/env python

"""

"""

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import numpy as np
import statistics
from numpy import linalg as la
from icecube.phys_services import I3Calculator
import operator
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.OpticalParameters import c, n, ngroup, theta_c


class ClusterPulseCleaning(icetray.I3ConditionalModule):
    """
    Causal Pulse Cleaning looks at the brightest DOM and assigns a time and location such that all other
    hits must be within a loose causality window. This Window is based on the travel of the particle at
    c + at most 50m of photon travel time at velocity c/n.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        # self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("inputseries", "GCD to be simulated", "MCPESeriesMap")
        self.AddParameter("output", "GCD to be simulated", "clusterpulses")
        self.AddParameter("cutsigma", "method to select cluster", 1.0)
        self.AddParameter("cutaveragesigma", "", 2.0)
        self.AddOutBox("OutBox")

    def Configure(self):
        # self.gcdFile = self.GetParameter("GCDFile")
        self.input = self.GetParameter("inputseries")
        self.output = self.GetParameter("output")
        self.cutsigma = self.GetParameter("cutsigma")
        self.cutaveragesigma = self.GetParameter("cutaveragesigma")

    # Compute base time and position that pulses need to be causal
    # with if no other pulses are causal with it.

    def Geometry(self, frame):
        self.domsUsed = frame["I3Geometry"].omgeo
        self.PushFrame(frame)

    def Physics(self, frame):
        mcpeMap = frame[self.input]
        clusterMCPEMap = dataclasses.I3RecoPulseSeriesMap()
        if len(mcpeMap.keys()) < 1:
            # self.PushFrame(frame)
            return
        removed = 10

        while removed > 0:
            distances = []
            distance_dict = {}
            n_under_200_dict = {}

            nopmtdomkeys = set()
            for omkey in mcpeMap.keys():
                nopmtdomkeys.add(NoPMTKey(omkey))

            for i, omkey1 in enumerate(nopmtdomkeys):
                for j, omkey2 in enumerate(nopmtdomkeys):
                    if j <= i:
                        continue
                    if NoPMTKey(omkey1) == NoPMTKey(omkey2):
                        continue
                    pos1 = self.domsUsed[AddPMTKey(omkey1, 1)].position
                    pos2 = self.domsUsed[AddPMTKey(omkey2, 1)].position
                    nopmtkey1 = NoPMTKey(omkey1)
                    nopmtkey2 = NoPMTKey(omkey2)

                    dom_to_dom = np.sqrt(
                        (pos1.x - pos2.x) ** 2.0
                        + (pos1.y - pos2.y) ** 2.0
                        + (pos1.z - pos2.z) ** 2.0
                    )
                    distances.append(dom_to_dom)

                    if dom_to_dom < 200.0:
                        if nopmtkey1 not in n_under_200_dict.keys():
                            n_under_200_dict[nopmtkey1] = 0
                        if nopmtkey2 not in n_under_200_dict.keys():
                            n_under_200_dict[nopmtkey2] = 0
                        n_under_200_dict[nopmtkey1] += 1
                        n_under_200_dict[nopmtkey2] += 1

                    if nopmtkey1 not in distance_dict.keys():
                        distance_dict[nopmtkey1] = dom_to_dom
                    else:
                        distance_dict[nopmtkey1] = min(
                            distance_dict[nopmtkey1], dom_to_dom
                        )
                    if nopmtkey2 not in distance_dict.keys():
                        distance_dict[nopmtkey2] = dom_to_dom
                    else:
                        distance_dict[nopmtkey2] = min(
                            distance_dict[nopmtkey2], dom_to_dom
                        )

            if len(distances) < 3:
                return

            mean = statistics.mean(distances)
            sigma = statistics.stdev(distances)

            removed = 0
            for omkey in mcpeMap.keys():
                # if distance_dict[NoPMTKey(omkey)] > mean + self.cutsigma*sigma :
                #    continue
                # print(distance_dict[NoPMTKey(omkey)])
                nopmtkey = NoPMTKey(omkey)
                if distance_dict[nopmtkey] > 100.0:
                    removed += 1
                    continue
                if n_under_200_dict[nopmtkey] < 3:
                    removed += 1
                    continue
                clusterMCPEMap[omkey] = mcpeMap[omkey]
            mcpeMap = clusterMCPEMap
            clusterMCPEMap = dataclasses.I3RecoPulseSeriesMap()

        frame[self.output] = mcpeMap
        self.PushFrame(frame)
