from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars
from SimFuncs import generateNoiseMCPEList

class addNoise(icetray.I3ConditionalModule):
    """
    Add Noise to the simulation as PEs
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile",
                          "GCD to be simulated",
                          '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
        self.AddParameter("PhysicsMCPETreeName",
                          "Name of the Physics I3MCTree name",
                          "MCPESeriesMap")
        self.AddParameter("NoiseMCPEMCTreeName",
                          "Name of the noise I3MCTree name",
                          "NoiseSeriesMap")
        self.AddParameter("STRAWTimestamps",
                          "Noise Hits for mDOMs to be injected")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.gcdFile = self.GetParameter("GCDFile")
        self.mcpeSeriesName = self.GetParameter("PhysicsMCPETreeName")
        self.noiseSeriesName = self.GetParameter("NoiseMCPEMCTreeName")
        #self.mDOM_Noise = self.GetParameter("mDOM_noise")
        self.timestamps = self.GetParameter("STRAWTimestamps")


    def DAQ(self, frame):

        gcd_file = dataio.I3File(self.gcdFile)
        cframe = gcd_file.pop_frame()
        geometry = cframe["I3Geometry"]
        geoMap = geometry.omgeo

        mcpeMap = frame[self.mcpeSeriesName]
        mcpeOMKeys = mcpeMap.keys()
        mcpeTime = ([])
        for mcpeList in mcpeMap.values():
            timeList = [mcpe.time for mcpe in mcpeList]
            mcpeTime = np.append(mcpeTime, timeList)
        medVal = np.mean(mcpeTime)
        #print 'mcpeTime', mcpeTime
        #print 'mean', medVal
        print 'median calculated'

        mcpeNoiseMap = simclasses.I3MCPESeriesMap()
        for modkey in geoMap.keys():
            #print 'generating noise in DOM'
            noiseList = generateNoiseMCPEList(modkey, medVal, self.timestamps)
            omkey = OMKey(modkey.string, modkey.om, 0)
            print omkey
            mcpeNoiseMap[omkey] = noiseList

        frame[self.noiseSeriesName] = mcpeNoiseMap

        self.PushFrame(frame)
