from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
import argparse
from SimFuncs import generateMCPEList

class makeHits(icetray.I3ConditionalModule):
    """
    Generate hits from I3Photons for P-ONE
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("WavelengthAcceptance",
                          "Path to the wavelength acceptance of the DOM",
                          "/home/users/akatil/P-ONE/git/PONE_NuTau/DOM/MDOM/DOMEfficiency.dat")
        self.AddParameter("AngularAcceptance",
                          "Input the flat angular acceptance for the mDOM",
                          0.811397114255)
        self.AddParameter("GCDFile",
                          "GCD to be simulated",
                          '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
        self.AddParameter("PropagatedPhotons",
                          "I3MCTree containing the list of propagated I3Photons",
                          "I3Photons")
        self.AddParameter("OutputMCTreeName",
                          "Name of the out put I3MCTree series name",
                          "MCPESeriesMap")
        self.AddParameter("DOMOversizeFactor",
                          "Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the CLSim run.\n"
	                      "The photon arrival times will be corrected. In practice this means your large spherical DOMs will\n"
	                      "become ellipsoids.",
                          1.)
        self.AddParameter("DefaultRelativeDOMEfficiency",
	                      "Default relative efficiency.",
	                      1.0);
        self.AddOutBox("OutBox")

    def Configure(self):

        self.wavelengthAcceptance = self.GetParameter("WavelengthAcceptance")
        self.angularAcceptance = self.GetParameter("AngularAcceptance")
        self.gcdFile = self.GetParameter("GCDFile")
        self.propagatedPhotons = self.GetParameter("PropagatedPhotons")
        self.outputMCTreeName = self.GetParameter("OutputMCTreeName")
        self.domOversizeFactor = self.GetParameter("DOMOversizeFactor")
        self.relativeEfficieny = self.GetParameter("DefaultRelativeDOMEfficiency")


    def DAQ(self, frame):

        gcd_file = dataio.I3File(self.gcdFile)
        cframe = gcd_file.pop_frame()
        geometry = cframe["I3Geometry"]
        geoMap = geometry.omgeo
        #calibration = cframe["I3Calibration"]
        #calMap = calibration.dom_cal

        DOMRadius = 0.16510*I3Units.m
        domAcceptance = clsim.GetIceCubeDOMAcceptance(
            	                                domRadius = DOMRadius*self.domOversizeFactor,
            	                                efficiency=self.relativeEfficieny)


        photonDOMMap = frame[self.propagatedPhotons]
        mcpeMap = simclasses.I3MCPESeriesMap()
        succPhotonMap = simclasses.I3CompressedPhotonSeriesMap()
        for modkey in photonDOMMap.keys():
            mcpeList, photonList = generateMCPEList(geoMap, photonDOMMap[modkey], modkey, self.angularAcceptance, domAcceptance, self.relativeEfficieny)
            omkey = OMKey(modkey.string, modkey.om, 0)
            if len(mcpeList) > 0:
                mcpeMap[omkey] = mcpeList
                succPhotonMap[modkey] = photonList

        if len(mcpeMap) > 0:
            frame[self.outputMCTreeName] = mcpeMap

            frame["SuccPhotonMap"] = succPhotonMap
            #frame.Delete("I3Photons")

            self.PushFrame(frame)
