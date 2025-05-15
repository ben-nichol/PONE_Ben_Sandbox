from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, DOMProperties, Geant4PMTAcceptance

# from DOMUtility import NoPMTKey, AddPMTKey, DOMProperties

import sys
sys.path.append('/home/jakubs/projects/def-mdanning/jakubs/k40/utils')
from POMModel import POM


class SimpleDOMSimulation(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("inputmap", "Name of the I3Photons from clsim", "I3Photons")
        self.AddParameter(
            "outputmap",
            "Name of the output I3RecoPulseSeriesMap",
            "I3RecoPulseSeriesMap",
        )
        self.AddParameter(
            "outputmap_mcpe", "Name of the output I3MCPESeriesMap", "I3MCPESeriesMap"
        )
        self.AddParameter("add_noise", "Should random noise be added?", True)
        self.AddParameter("PMT_tts", "Transit time spread of PMT", 3.0 * I3Units.ns)
        self.AddParameter("PMT_ts", "Transit time of PMT", 25.0 * I3Units.ns)
        self.AddParameter("chargesigma", "Sigma of charge distribution", 0.3)
        self.AddParameter("chargemean", "Mean of Charge distribution ", 1.0)
        self.AddParameter("DNprob", "Dark Noise rate (pulses per ns)", 0.000001)
        self.AddParameter("APprob", "Total AP probability", 0.06)
        self.AddParameter(
            "APmeantime_1", "Mean of early time AP distribution", 2000.0 * I3Units.ns
        )
        self.AddParameter(
            "APtimesigma_1", "Sigma of early time AP distribution", 1000.0 * I3Units.ns
        )
        self.AddParameter(
            "APmeantime_2", " mean of late time AP distribution", 8000.0 * I3Units.ns
        )
        self.AddParameter(
            "APtimesigma_2", " sigma of late time AP distribution", 2000.0 * I3Units.ns
        )
        self.AddParameter(
            "APComponetRatio", "fraction of AP in early time component", 0.3
        )
        self.AddParameter("LPprob", "Probability for a late pulse", 0.01)
        self.AddParameter("minTsep", " minimum time for a separated pulse.", 3.0)
        self.AddParameter("PEthreshold", " Pulse charge threshold", 0.25)
        self.AddParameter("PEsaturation", " Saturation threshold for PMT", 100.0)
        self.AddParameter("RandomService", "Random Service")
        self.AddParameter("SplitDoms", "", True)
        self.AddParameter("DOMAcceptanceFile", "", "")
        self.AddParameter("PMTQEFile", "", "")
        self.AddParameter("DropStrings", "", [])
        self.AddParameter("NoisePulseSeries", "", [])
        self.AddParameter("NoPureNoiseEvents", "", True)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.inputmap = self.GetParameter("inputmap")
        self.outputmap = self.GetParameter("outputmap")
        self.outputmap_mcpe = self.GetParameter("outputmap_mcpe")
        self.add_noise = self.GetParameter("add_noise")
        self.PMT_tts = self.GetParameter("PMT_tts")
        self.PMT_ts = self.GetParameter("PMT_ts")
        self.chargesigma = self.GetParameter("chargesigma")
        self.chargemean = self.GetParameter("chargemean")
        self.DNprob = self.GetParameter("DNprob")
        self.APprob = self.GetParameter("APprob")
        self.APmeantime_1 = self.GetParameter("APmeantime_1")
        self.APtimesigma_1 = self.GetParameter("APtimesigma_1")
        self.APmeantime_2 = self.GetParameter("APmeantime_2")
        self.APtimesigma_2 = self.GetParameter("APtimesigma_2")
        self.APComponetRatio = self.GetParameter("APComponetRatio")
        self.LPprob = self.GetParameter("LPprob")
        self.minTsep = self.GetParameter("minTsep")
        self.PEthreshold = self.GetParameter("PEthreshold")
        self.PEsaturation = self.GetParameter("PEsaturation")
        self.randomService = self.GetParameter("RandomService")
        self.splitDOMs = self.GetParameter("SplitDoms")
        self.dropstrings = self.GetParameter("DropStrings")
        self.noisePulseSeries = self.GetParameter("NoisePulseSeries")
        self.noPureNoiseEvents = self.GetParameter("NoPureNoiseEvents")
        if self.GetParameter("PMTQEFile") != "":
            GetPMTQETable(self.photonweights, self.GetParameter("PMTQEFile"))

        kwargs_DOMProperties = {}
        if self.GetParameter("DOMAcceptanceFile") != "":
            kwargs_DOMProperties["PMTAcceptanceFile"] = self.GetParameter(
                "DOMAcceptanceFile"
            )
        if self.GetParameter("PMTQEFile") != "":
            kwargs_DOMProperties["PMTQEFile"] = self.GetParameter("PMTQEFile")

        # load the DOM properties from their respective configuration files
        self.dom_properties = DOMProperties(**kwargs_DOMProperties)
        self.WroteActiveDOMsToSimFrame = False

    #########################################################################

    # split Photon hits into PMTs on the DOMs.
    def SplitPMTs(self, photonmap, dropstrings=[]):
        newphotonmap = {}
        pulseseries = dataclasses.I3RecoPulseSeries()
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        # make new map with individual PMTs

        for omkey in photonmap.keys():
            if omkey.string in dropstrings:
                continue
            newomkey = NoPMTKey(omkey)
            newphotonmap[newomkey] = []
            for pulse in photonmap[omkey]:
                pmtid = self.GetPMT(
                    photonDir=[pulse.dir.x, pulse.dir.y, pulse.dir.z],
                    wl=pulse.wavelength / I3Units.nanometer,
                    weight=pulse.weight,
                )

                if pmtid < 0:
                    continue
                newphotonmap[newomkey].append((pulse.time, pmtid))
                pmtkey = AddPMTKey(newomkey, pmtid)
                if pmtkey not in outputpulsemap.keys():
                    outputpulsemap[pmtkey] = dataclasses.I3RecoPulseSeries()
                rpulse = dataclasses.I3RecoPulse()
                rpulse.time = pulse.time
                rpulse.charge = 1.0
                outputpulsemap[pmtkey].append(rpulse)

        # print("split")
        # print(newphotonmap)
        return newphotonmap, outputpulsemap

    """!
    GetPMT(photonDir,wl,random)
    Input: 
        photonDir = list of x,y,z or theta,phi
        wl = wavelength (m or nm)
        random = random number between 0 and 1
    Operation:
        Randomly assigns photon hits to PMTs based on the photon direction and PMT acceptance.
        This function also applies the QE shape. This function assumes that CLSim has already
        scaled down the photon production so that the maxacceptance after CLSim is scaled to 1.0.
    """

    def GetPMT(self, photonDir, wl, weight):
        random = self.randomService.uniform(0.0, 1.0)

        theta = 0.0
        phi = 0.0
        if len(photonDir) < 3:
            theta = photonDir[0]
            phi = photonDir[1]
        else:
            theta = np.arccos(photonDir[2])
            phi = np.arctan2(photonDir[1], photonDir[0])

        while phi < 0.0:
            phi += 2.0 * np.pi

        thetaBin = max(0, min(178, int(180.0 * theta / np.pi)))
        phiBin = max(0, min(358, int(180.0 * phi / np.pi)))
        pmtprobs = []
        for i in range(len(self.dom_properties.PMTacceptance)):
            pmtprobs.append(
                self.dom_properties.PMTacceptance[i][thetaBin][phiBin]
                / self.dom_properties.maxAngularAcceptance
            )

        totalprob = sum(pmtprobs)
        if random > totalprob:
            # print("angular kill photon")
            return -1

        random = self.randomService.uniform(0.0, 1.0)
        i = 0
        sumprob = pmtprobs[0]
        while random > sumprob / totalprob and i < len(pmtprobs) - 1:
            i += 1
            sumprob += pmtprobs[i]

        qe = self.dom_properties.GetPMTQEnm(wl)
        # note: the probaility should only be weight * qe. It also includes maxAngularAcceptance
        #       here, since we already took this factor into account in `GetPMT` (which also rejects
        #       hits.)
        probability = weight * qe * self.dom_properties.maxAngularAcceptance

        # print("probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format( probability, weight, qe, self.dom_properties.maxAngularAcceptance, wl ))

        if probability > 1.0:
            print(
                "ERROR: probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format(
                    probability,
                    weight,
                    qe,
                    self.dom_properties.maxAngularAcceptance,
                    wl,
                )
            )
            raise RuntimeError(
                "The combined detection probability should never be > 1. You need to re-generate your I3Photons with a higher overall bias weight (setting `WavelengthAcceptance`)"
            )

        ### reject photon according to `probability`
        random = self.randomService.uniform(0.0, 1.0)
        if random > probability:
            # print("kill photon")
            return -1

        # print("pmt = "+str(i+1))

        return i + 1

    # Get the min and max pulse times to set time window for dark hits.
    def GetMaxMinTimes(self, mcpemap):
        max_pt = -999999999.0
        min_pt = 999999999.0
        # print(mcpemap)
        for omkey in mcpemap.keys():
            for pulse in mcpemap[omkey]:
                if max_pt < pulse[0]:
                    max_pt = pulse[0]
                if min_pt > pulse[0]:
                    min_pt = pulse[0]

        max_pt += 10000.0
        min_pt -= 2000.0

        return max_pt, min_pt

    # Add Dark hits across all DOMs (addes I3MCPEs)
    def AddDarkHits(self, max_pt, min_pt):
        newmcpemap = {}
        # print("min = "+str(min_pt)+" max = "+str(max_pt))
        npmts = int(self.dom_properties.GetNPMTs())
        for omkey in self.domkeys:
            # print(omkey)
            dt = self.randomService.exp(1.0 / (self.DNprob * npmts)) + min_pt
            ipmt = 1 + self.randomService.integer(npmts)
            first = True
            while dt < max_pt:
                if first:
                    newmcpemap[omkey] = []
                    first = False
                # print(dt)
                newmcpemap[omkey].append((dt, ipmt))
                dt += self.randomService.exp(1.0 / (self.DNprob * npmts))
                ipmt = 1 + self.randomService.integer(npmts)

        return newmcpemap

    def AddEnvironmentNoise(self, darkhits, noisepulses):
        newmcpemap = {}
        mergedmap = {}
        for pulseseries in noisepulses:
            if type(pulseseries) == type(dataclasses.I3RecoPulseSeries()):
                for dom in pulseseries.keys():
                    nopmtkey = NoPMTKey(dom)
                    if nopmtkey not in newmcpemap.keys():
                        newmcpemap[nopmtkey] = list()
                    newmcpemap[nopmtkey].extend(
                        [
                            (pulseseries[dom][i].time, dom.pmt)
                            for i in range(len(pulseseries[dom]))
                        ]
                    )
            elif type(pulseseries) == type(simclasses, I3PhotonSeries()):
                for dom in noisepulses.keys():
                    newmcpemap[nopmtkey(dom)] = list()
                    for pulse in noisepulses[dom]:
                        pmtid = self.GetPMT(
                            photonDir=[pulse.dir.x, pulse.dir.y, pulse.dir.z],
                            wl=pulse.wavelength / I3Units.nanometer,
                            weight=pulse.weight,
                        )
                        newmcpemap[nopmtkey].append((pulse.time, pmtid))
            else:
                print("Invalid noise pulse series")
        for dom in newmcpemap.keys():
            newmcpemap[dom].sort(key=lambda x: x[0])

        for dom in newmcpemap.keys():
            mergedmap[dom] = list()
            if dom in darkhits.keys():
                i = 0
                j = 0
                while i < len(darkhits[dom]) and j < len(newmcpemap[dom]):
                    if darkhits[dom].time < sortedpulses[j].time:
                        mergedmap[dom].append(darkhits[dom])
                        i += 1
                    else:
                        mergedmap[dom].append(sortedpulses[j])
                        j += 1
                if i < len(darkhits[dom]):
                    mergedmap[dom].append(darkhits[dom])
                    i += 1
                if j < len(sortedpulses):
                    mergedmap[dom].append(sortedpulses[j])
                    j += 1
            else:
                mergedmap[dom] = newmcpemap[dom]

        for dom in darkhits.keys():
            if dom not in mergedmap.keys():
                mergedmap[dom] = darkhits[dom]
        return mergedmap

    def CombineOrderedLists(self, list1, list2):
        pulsetimelist = []
        j = 0
        i = 0
        while j < len(list1) and i < len(list2):
            if list1[j][0] <= list2[i][0]:
                pulsetimelist.append(list1[j])
                j += 1
            else:
                pulsetimelist.append(list2[i])
                i += 1
        while j < len(list1):
            pulsetimelist.append(list1[j])
            j += 1
        while i < len(list2):
            pulsetimelist.append(list2[i])
            i += 1
        return pulsetimelist

    def ApplyPMTCharacteristics(self, mcpemap):
        pulsetimelist = []
        pulseseries = dataclasses.I3RecoPulseSeries()
        dompulseseries = dataclasses.I3RecoPulseSeries()
        pulseoutoftimelist = []
        for pe in mcpemap:
            time = self.randomService.gaus(pe[0], self.PMT_tts)
            if self.randomService.uniform(0.0, 1.0) < self.LPprob:
                time += self.randomService.gaus(
                    self.PMT_ts * 2.0, np.sqrt(2.0) * self.PMT_tts
                )
            if len(pulsetimelist) < 1 or time > pulsetimelist[-1][0]:
                pulsetimelist.append((time, pe[1]))
            else:
                pulseoutoftimelist.append((time, pe[1]))
        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(
                pulsetimelist.copy(), pulseoutoftimelist
            )

        for pe in pulsetimelist:
            pulseoutoftimelist = []
            if self.randomService.uniform(0.0, 1.0) < self.APprob:
                if self.randomService.uniform(0.0, 1.0) < self.APComponetRatio:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_1, self.APtimesigma_1
                    )
                else:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_2, self.APtimesigma_2
                    )
                pulseoutoftimelist.append((time, pe[1]))

        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(
                pulsetimelist.copy(), pulseoutoftimelist
            )

        return pulsetimelist

    def ApplyPMTCharacteristicsDarkNoise(self, pelist):
        pulseoutoftimelist = []
        pulsetimelist = []

        for pe in pelist:
            pulseoutoftimelist = []
            if self.randomService.uniform(0.0, 1.0) < self.APprob:
                if self.randomService.uniform(0.0, 1.0) < self.APComponetRatio:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_1, self.APtimesigma_1
                    )
                else:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_2, self.APtimesigma_2
                    )
                pulseoutoftimelist.append((time, pe[1]))

        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(pelist, pulseoutoftimelist)
        else:
            pulsetimelist = pelist

        return pulsetimelist

    def MakeRecoPulses(
        self, pulsetimelist, pulsechargelist, omkey, outputpulsemap, DOMpulsemap=None
    ):
        mingap = 4.0
        minindex = -1

        if len(pulsetimelist) > 100:
            leading = 0
            following = 1
            while following < len(pulsetimelist):
                if (
                    pulsetimelist[following][0] - pulsetimelist[leading][0]
                ) < 3.0 and pulsechargelist[leading] * pulsechargelist[following] > 0.0:
                    pulsechargelist[leading] += pulsechargelist[following]
                    pulsechargelist[following] = 0.0
                elif pulsechargelist[following] > 0.0:
                    leading = following
                following += 1
        else:
            # needs to be better
            for i in range(1, len(pulsetimelist)):
                if (
                    pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                ) < mingap and pulsechargelist[i] * pulsechargelist[i - 1] > 0.0:
                    mingap = pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                    minindex = i
            # If less than limit, combine pulses
            while mingap <= self.minTsep:
                if pulsechargelist[minindex] > pulsechargelist[minindex - 1]:
                    pulsechargelist[minindex] += pulsechargelist[minindex - 1]
                    pulsechargelist[minindex - 1] = 0.0
                else:
                    pulsechargelist[minindex - 1] += pulsechargelist[minindex]
                    pulsechargelist[minindex] = 0.0
                mingap = self.minTsep + 1.0
                minindex = -1
                # reestablish new min gap
                for i in range(1, len(pulsetimelist)):
                    if (
                        pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                    ) < mingap and pulsechargelist[i] * pulsechargelist[i - 1] > 0.0:
                        mingap = pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                        minindex = i

        pmtinlist = []
        for i in range(len(pulsetimelist)):
            if pulsechargelist[-1 - i] < self.PEthreshold:
                continue
            if pulsetimelist[i][1] not in pmtinlist:
                pmtinlist.append(pulsetimelist[i][1])
        pmtinlist.sort()

        for i in range(len(pmtinlist)):
            newomkey = AddPMTKey(omkey, pmtinlist[i])
            outputpulsemap[newomkey] = dataclasses.I3RecoPulseSeries()

        for i in range(len(pulsetimelist)):
            # remove pulses with too low charge.
            if pulsechargelist[-1 - i] < self.PEthreshold:
                continue
            rpulse = dataclasses.I3RecoPulse()
            rpulse.time = pulsetimelist[i][0]
            # saturate pulses with too much charge.
            if pulsechargelist[-1 - i] > self.PEsaturation:
                rpulse.charge = self.PEsaturation
                rpulse.charge += (pulsechargelist[-1 - i] - self.PEsaturation) * (
                    self.PEsaturation / pulsechargelist[-1 - i]
                )
            else:
                rpulse.charge = pulsechargelist[-1 - i]
            rpulse.width = pulsetimelist[i][1]
            if not (DOMpulsemap is None):
                if omkey not in DOMpulsemap.keys():
                    DOMpulsemap[omkey] = dataclasses.I3RecoPulseSeries()
                DOMpulsemap[omkey].append(rpulse)
            newomkey = AddPMTKey(omkey, pulsetimelist[i][1])
            outputpulsemap[newomkey].append(rpulse)

    # Apply the responce of the PMT, this includes pulse combining for pulses too close together.
    def ApplyPMTResponse(self, mcpemap, darkhits):
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        outputpulsemap_nonoise = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsemap = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsetimelist = {}

        self.npmt = 0

        for istring in range(1, self.nstring + 1):
            for iom in range(1, self.nom + 1):
                omkey = OMKey(istring, iom, 0)
                if omkey in mcpemap.keys():
                    RealPulses = self.ApplyPMTCharacteristics(mcpemap[omkey])
                    # print("real pulse = "+str(len(mcpemap[omkey])))
                    pulsetimelist = []
                    realpulsetimelist = []
                    isreal = []
                    DarkPulses = []
                    if omkey in darkhits.keys():
                        DarkPulses = self.ApplyPMTCharacteristicsDarkNoise(
                            darkhits[omkey]
                        )
                        # print("dark")
                        # print(len(DarkPulses))
                        # print(len(darkhits[omkey]))
                        j = 0
                        i = 0
                        while j < len(RealPulses) and i < len(DarkPulses):
                            if RealPulses[j][0] <= DarkPulses[i][0]:
                                pulsetimelist.append(RealPulses[j])
                                isreal.append(1.0)
                                j += 1
                            else:
                                pulsetimelist.append(DarkPulses[i])
                                isreal.append(0.0)
                                i += 1
                        while j < len(RealPulses):
                            pulsetimelist.append(RealPulses[j])
                            isreal.append(1.0)
                            j += 1
                        while i < len(DarkPulses):
                            pulsetimelist.append(DarkPulses[i])
                            isreal.append(0.0)
                            i += 1
                    else:
                        pulsetimelist = RealPulses
                        for i in range(len(pulsetimelist)):
                            isreal.append(1.0)

                    pulsechargelist = []
                    realpulsechargelist = []

                    for i in range(len(pulsetimelist)):
                        charge = self.randomService.gaus(
                            self.chargemean, self.chargesigma
                        )
                        pulsechargelist.append(charge)
                        if isreal[i] > 0.0:
                            realpulsechargelist.append(charge)

                    # print(len(RealPulses))
                    # print(len(realpulsechargelist))
                    # print(len(pulsetimelist))
                    # print(len(pulsechargelist))
                    # print("MakeRecos")

                    if len(DarkPulses) > 0 and len(RealPulses) > 0:
                        self.MakeRecoPulses(
                            RealPulses,
                            realpulsechargelist,
                            omkey,
                            outputpulsemap_nonoise,
                        )
                        self.MakeRecoPulses(
                            pulsetimelist,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )
                    elif len(RealPulses) > 0:
                        self.MakeRecoPulses(
                            RealPulses,
                            realpulsechargelist,
                            omkey,
                            outputpulsemap_nonoise,
                            DOMpulsemap,
                        )
                        for ipmt in range(1, self.dom_properties.GetNPMTs() + 1):
                            newomkey = AddPMTKey(omkey, ipmt)
                            if newomkey in outputpulsemap_nonoise.keys():
                                outputpulsemap[newomkey] = outputpulsemap_nonoise[
                                    newomkey
                                ]
                    elif len(DarkPulses) > 0:
                        self.MakeRecoPulses(
                            pulsetimelist,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )
                    else:
                        continue

                elif omkey in darkhits.keys():
                    if omkey in mcpemap.keys():
                        continue
                    if len(darkhits[omkey]) > 0:
                        DarkPulses = self.ApplyPMTCharacteristicsDarkNoise(
                            darkhits[omkey]
                        )
                        pulsechargelist = []
                        for i in range(len(DarkPulses)):
                            pulsechargelist.append(
                                self.randomService.gaus(
                                    self.chargemean, self.chargesigma
                                )
                            )

                        self.MakeRecoPulses(
                            DarkPulses,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )

                # if omkey in darkhits.keys() or omkey in mcpemap.keys() :
                #    for ipmt in range(self.dom_properties.GetNPMTs()) :
                #        newomkey = AddPMTKey(omkey,ipmt)
                #        print(newomkey)
                #        print("all")
                #        if  newomkey in outputpulsemap.keys() :
                #            print(outputpulsemap[newomkey])
                #        else:
                #            print("no pulses")
                #        print("no noise")
                #        if newomkey in outputpulsemap_nonoise.keys() :
                #            print(outputpulsemap_nonoise[newomkey])
                #        else :
                #            print("no pulses")

        return outputpulsemap, outputpulsemap_nonoise, DOMpulsemap

    def Geometry(self, frame):
        if len(self.dropstrings) > 0:
            self.domsUsed = dataclasses.I3OMGeoMap()
            for omkey in frame["I3Geometry"].omgeo.keys():
                if omkey.string in self.dropstrings:
                    continue
                self.domsUsed[omkey] = frame["I3Geometry"].omgeo[omkey]
        else:
            self.domsUsed = frame["I3Geometry"].omgeo

        domkeylist = []
        self.nstring = 0
        self.nom = 0
        self.npmt = 0
        for omkey in self.domsUsed.keys():
            self.nstring = max(self.nstring, omkey.string)
            self.nom = max(self.nom, omkey.om)
            self.npmt = max(self.npmt, omkey.pmt)
            domkeylist.append(NoPMTKey(omkey))
        self.domkeys = set(domkeylist)
        self.domkeys = sorted(self.domkeys, key=lambda x: (x.string, x.om, x.pmt))
        self.nstring += 1
        self.nom += 1
        self.npmt += 1
        self.PushFrame(frame)

    #########################################################################

    def Simulation(self, frame):
        if len(self.dropstrings) > 0:
            frame["SimulatedDOMs"] = self.domsUsed

        self.WroteActiveDOMsToSimFrame = True
        self.PushFrame(frame)

    def DAQ(self, frame):
        if not self.WroteActiveDOMsToSimFrame:
            simframe = icetray.I3Frame("S")
            self.Simulation(simframe)

        photonmap = frame[self.inputmap]
        noisepulses = [
            frame[noisepulsename]
            for noisepulsename in self.noisePulseSeries
            if frame.Has(noisepulsename)
        ]
        lennoisepulses = sum([len(noisepulses[i]) for i in range(len(noisepulses))])
        # print("ndoms = "+str(len(photonmap)))
        # split photons from DOMs to PMTs-on-DOMs
        # and apply weights
        if (len(photonmap) < 1) and self.noPureNoiseEvents:
            return
        if len(photonmap) < 1 and lennoisepulses < 1:
            return

        #print("Splitting PMTs")
        photonmap_on_pmts, frame[self.inputmap + "_pmtsplit"] = self.SplitPMTs(
            photonmap, self.dropstrings
        )
        # print(" pmt split pulses = " + str(len(photonmap_on_pmts)))
        # find the minimum and maximum time (with padding)
        max_pt, min_pt = self.GetMaxMinTimes(photonmap_on_pmts)
        # print("get min max")
        darkhits = {}
        if self.add_noise:
            # add dark noise
            # print("add noise")
            darkhits = self.AddDarkHits(max_pt, min_pt)
            darkhits = self.AddEnvironmentNoise(darkhits, noisepulses)
        # print("ndarkhitpmts = " + str(len(darkhits)))
        # Save the MCPEs to the frame
        # frame[self.outputmap_mcpe] = mcpemap
        # print("apply pmt response")
        (
            frame[self.outputmap],
            frame[self.outputmap + "_nonoise"],
            frame["triggerpulsemap"],
        ) = self.ApplyPMTResponse(photonmap_on_pmts, darkhits)
        # print(frame['triggerpulsemap'])
        # print("final")
        # print(frame[self.outputmap])

        self.PushFrame(frame)


class SimpleDOMSimulationNew(icetray.I3ConditionalModule):
    """
    Drop in replacement of SimpleDOMSimulation that uses Geant4PMTAcceptance.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("inputmap", "Name of the I3Photons from clsim", "I3Photons")
        self.AddParameter(
            "outputmap",
            "Name of the output I3RecoPulseSeriesMap",
            "I3RecoPulseSeriesMap",
        )
        self.AddParameter(
            "outputmap_mcpe", "Name of the output I3MCPESeriesMap", "I3MCPESeriesMap"
        )
        self.AddParameter("add_noise", "Should random noise be added?", True)
        self.AddParameter("PMT_tts", "Transit time spread of PMT", 3.0 * I3Units.ns)
        self.AddParameter("PMT_ts", "Transit time of PMT", 25.0 * I3Units.ns)
        self.AddParameter("chargesigma", "Sigma of charge distribution", 0.3)
        self.AddParameter("chargemean", "Mean of Charge distribution ", 1.0)
        self.AddParameter("DNprob", "Dark Noise rate (pulses per ns)", 0.000001)
        self.AddParameter("APprob", "Total AP probability", 0.06)
        self.AddParameter(
            "APmeantime_1", "Mean of early time AP distribution", 2000.0 * I3Units.ns
        )
        self.AddParameter(
            "APtimesigma_1", "Sigma of early time AP distribution", 1000.0 * I3Units.ns
        )
        self.AddParameter(
            "APmeantime_2", " mean of late time AP distribution", 8000.0 * I3Units.ns
        )
        self.AddParameter(
            "APtimesigma_2", " sigma of late time AP distribution", 2000.0 * I3Units.ns
        )
        self.AddParameter(
            "APComponetRatio", "fraction of AP in early time component", 0.3
        )
        self.AddParameter("LPprob", "Probability for a late pulse", 0.01)
        self.AddParameter("minTsep", " minimum time for a separated pulse.", 3.0)
        self.AddParameter("PEthreshold", " Pulse charge threshold", 0.25)
        self.AddParameter("PEsaturation", " Saturation threshold for PMT", 100.0)
        self.AddParameter("RandomService", "Random Service")
        self.AddParameter("SplitDoms", "", True)
        self.AddParameter("DOMAcceptanceFile", "", "")
        self.AddParameter("PMTQEFile", "", "")
        self.AddParameter("DropStrings", "", [])
        self.AddParameter("NoisePulseSeries", "", [])
        self.AddParameter("NoPureNoiseEvents", "", True)

        self.pmt_acc = Geant4PMTAcceptance()

    def Configure(self):
        self.inputmap = self.GetParameter("inputmap")
        self.outputmap = self.GetParameter("outputmap")
        self.outputmap_mcpe = self.GetParameter("outputmap_mcpe")
        self.add_noise = self.GetParameter("add_noise")
        self.PMT_tts = self.GetParameter("PMT_tts")
        self.PMT_ts = self.GetParameter("PMT_ts")
        self.chargesigma = self.GetParameter("chargesigma")
        self.chargemean = self.GetParameter("chargemean")
        self.DNprob = self.GetParameter("DNprob")
        self.APprob = self.GetParameter("APprob")
        self.APmeantime_1 = self.GetParameter("APmeantime_1")
        self.APtimesigma_1 = self.GetParameter("APtimesigma_1")
        self.APmeantime_2 = self.GetParameter("APmeantime_2")
        self.APtimesigma_2 = self.GetParameter("APtimesigma_2")
        self.APComponetRatio = self.GetParameter("APComponetRatio")
        self.LPprob = self.GetParameter("LPprob")
        self.minTsep = self.GetParameter("minTsep")
        self.PEthreshold = self.GetParameter("PEthreshold")
        self.PEsaturation = self.GetParameter("PEsaturation")
        self.randomService = self.GetParameter("RandomService")
        self.splitDOMs = self.GetParameter("SplitDoms")
        self.dropstrings = self.GetParameter("DropStrings")
        self.noisePulseSeries = self.GetParameter("NoisePulseSeries")
        self.noPureNoiseEvents = self.GetParameter("NoPureNoiseEvents")
        if self.GetParameter("PMTQEFile") != "":
            GetPMTQETable(self.photonweights, self.GetParameter("PMTQEFile"))

        kwargs_DOMProperties = {}
        if self.GetParameter("DOMAcceptanceFile") != "":
            kwargs_DOMProperties["PMTAcceptanceFile"] = self.GetParameter(
                "DOMAcceptanceFile"
            )
        if self.GetParameter("PMTQEFile") != "":
            kwargs_DOMProperties["PMTQEFile"] = self.GetParameter("PMTQEFile")

        # load the DOM properties from their respective configuration files
        self.dom_properties = DOMProperties(**kwargs_DOMProperties)
        self.WroteActiveDOMsToSimFrame = False

    #########################################################################

    # split Photon hits into PMTs on the DOMs.
    def SplitPMTs(self, photonmap, dropstrings=[]):
        newphotonmap = {}
        pulseseries = dataclasses.I3RecoPulseSeries()
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        # make new map with individual PMTs

        for omkey in photonmap.keys():
            if omkey.string in dropstrings:
                continue
            newomkey = NoPMTKey(omkey)
            newphotonmap[newomkey] = []

            n_photons = len(photonmap[omkey])
            photon_positions = np.empty((n_photons, 3))
            photon_wavelengths = np.empty(n_photons)
            photon_weights = np.empty(n_photons)

            for i, photon in enumerate(photonmap[omkey]):

                photon_positions[i] = [photon.pos.x, photon.pos.y, photon.pos.z]
                photon_wavelengths[i] = photon.wavelength / I3Units.nanometer
                photon_weights[i] = photon.weight
            photon_positions /= self.modgeo[omkey].radius

            pmt_ids = self.pmt_acc.check_pmt_hit(
                photon_positions, photon_wavelengths, photon_weights
            )

            for photon, pmtid in zip(photonmap[omkey], pmt_ids):
                if pmtid == 0:
                    continue
                newphotonmap[newomkey].append((photon.time, pmtid))
                pmtkey = AddPMTKey(newomkey, int(pmtid))
                if pmtkey not in outputpulsemap.keys():
                    outputpulsemap[pmtkey] = dataclasses.I3RecoPulseSeries()
                rpulse = dataclasses.I3RecoPulse()
                rpulse.time = photon.time
                rpulse.charge = 1.0
                outputpulsemap[pmtkey].append(rpulse)

        # print("split")
        # print(newphotonmap)
        return newphotonmap, outputpulsemap

    # Get the min and max pulse times to set time window for dark hits.
    def GetMaxMinTimes(self, mcpemap):
        max_pt = -999999999.0
        min_pt = 999999999.0
        # print(mcpemap)
        for omkey in mcpemap.keys():
            for pulse in mcpemap[omkey]:
                if max_pt < pulse[0]:
                    max_pt = pulse[0]
                if min_pt > pulse[0]:
                    min_pt = pulse[0]

        max_pt += 10000.0
        min_pt -= 2000.0

        return max_pt, min_pt

    # Add Dark hits across all DOMs (addes I3MCPEs)
    def AddDarkHits(self, max_pt, min_pt):
        newmcpemap = {}
        # print("min = "+str(min_pt)+" max = "+str(max_pt))
        npmts = int(self.dom_properties.GetNPMTs())
        for omkey in self.domkeys:
            # print(omkey)
            dt = self.randomService.exp(1.0 / (self.DNprob * npmts)) + min_pt
            ipmt = 1 + self.randomService.integer(npmts)
            first = True
            while dt < max_pt:
                if first:
                    newmcpemap[omkey] = []
                    first = False
                # print(dt)
                newmcpemap[omkey].append((dt, ipmt))
                dt += self.randomService.exp(1.0 / (self.DNprob * npmts))
                ipmt = 1 + self.randomService.integer(npmts)

        return newmcpemap

    def AddEnvironmentNoise(self, darkhits, noisepulses):
        newmcpemap = {}
        mergedmap = {}
        for pulseseries in noisepulses:
            if type(pulseseries) == type(dataclasses.I3RecoPulseSeries()):
                for dom in pulseseries.keys():
                    nopmtkey = NoPMTKey(dom)
                    if nopmtkey not in newmcpemap.keys():
                        newmcpemap[nopmtkey] = list()
                    newmcpemap[nopmtkey].extend(
                        [
                            (pulseseries[dom][i].time, dom.pmt)
                            for i in range(len(pulseseries[dom]))
                        ]
                    )
            elif type(pulseseries) == type(simclasses, I3PhotonSeries()):
                for dom in noisepulses.keys():
                    newmcpemap[nopmtkey(dom)] = list()
                    for pulse in noisepulses[dom]:
                        pmtid = self.GetPMT(
                            photonDir=[pulse.dir.x, pulse.dir.y, pulse.dir.z],
                            wl=pulse.wavelength / I3Units.nanometer,
                            weight=pulse.weight,
                        )
                        newmcpemap[nopmtkey].append((pulse.time, pmtid))
            else:
                print("Invalid noise pulse series")
        for dom in newmcpemap.keys():
            newmcpemap[dom].sort(key=lambda x: x[0])

        for dom in newmcpemap.keys():
            mergedmap[dom] = list()
            if dom in darkhits.keys():
                i = 0
                j = 0
                while i < len(darkhits[dom]) and j < len(newmcpemap[dom]):
                    if darkhits[dom].time < sortedpulses[j].time:
                        mergedmap[dom].append(darkhits[dom])
                        i += 1
                    else:
                        mergedmap[dom].append(sortedpulses[j])
                        j += 1
                if i < len(darkhits[dom]):
                    mergedmap[dom].append(darkhits[dom])
                    i += 1
                if j < len(sortedpulses):
                    mergedmap[dom].append(sortedpulses[j])
                    j += 1
            else:
                mergedmap[dom] = newmcpemap[dom]

        for dom in darkhits.keys():
            if dom not in mergedmap.keys():
                mergedmap[dom] = darkhits[dom]
        return mergedmap

    def CombineOrderedLists(self, list1, list2):
        pulsetimelist = []
        j = 0
        i = 0
        while j < len(list1) and i < len(list2):
            if list1[j][0] <= list2[i][0]:
                pulsetimelist.append(list1[j])
                j += 1
            else:
                pulsetimelist.append(list2[i])
                i += 1
        while j < len(list1):
            pulsetimelist.append(list1[j])
            j += 1
        while i < len(list2):
            pulsetimelist.append(list2[i])
            i += 1
        return pulsetimelist

    def ApplyPMTCharacteristics(self, mcpemap):
        pulsetimelist = []
        pulseseries = dataclasses.I3RecoPulseSeries()
        dompulseseries = dataclasses.I3RecoPulseSeries()
        pulseoutoftimelist = []
        for pe in mcpemap:
            time = self.randomService.gaus(pe[0], self.PMT_tts)
            if self.randomService.uniform(0.0, 1.0) < self.LPprob:
                time += self.randomService.gaus(
                    self.PMT_ts * 2.0, np.sqrt(2.0) * self.PMT_tts
                )
            if len(pulsetimelist) < 1 or time > pulsetimelist[-1][0]:
                pulsetimelist.append((time, pe[1]))
            else:
                pulseoutoftimelist.append((time, pe[1]))
        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(
                pulsetimelist.copy(), pulseoutoftimelist
            )

        for pe in pulsetimelist:
            pulseoutoftimelist = []
            if self.randomService.uniform(0.0, 1.0) < self.APprob:
                if self.randomService.uniform(0.0, 1.0) < self.APComponetRatio:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_1, self.APtimesigma_1
                    )
                else:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_2, self.APtimesigma_2
                    )
                pulseoutoftimelist.append((time, pe[1]))

        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(
                pulsetimelist.copy(), pulseoutoftimelist
            )

        return pulsetimelist

    def ApplyPMTCharacteristicsDarkNoise(self, pelist):
        pulseoutoftimelist = []
        pulsetimelist = []

        for pe in pelist:
            pulseoutoftimelist = []
            if self.randomService.uniform(0.0, 1.0) < self.APprob:
                if self.randomService.uniform(0.0, 1.0) < self.APComponetRatio:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_1, self.APtimesigma_1
                    )
                else:
                    time = pe[0] + self.randomService.gaus(
                        self.APmeantime_2, self.APtimesigma_2
                    )
                pulseoutoftimelist.append((time, pe[1]))

        if len(pulseoutoftimelist) > 0:
            pulseoutoftimelist.sort(key=lambda x: x[0])
            pulsetimelist = self.CombineOrderedLists(pelist, pulseoutoftimelist)
        else:
            pulsetimelist = pelist

        return pulsetimelist

    def MakeRecoPulses(
        self, pulsetimelist, pulsechargelist, omkey, outputpulsemap, DOMpulsemap=None
    ):
        mingap = 4.0
        minindex = -1

        if len(pulsetimelist) > 100:
            leading = 0
            following = 1
            while following < len(pulsetimelist):
                if (
                    pulsetimelist[following][0] - pulsetimelist[leading][0]
                ) < 3.0 and pulsechargelist[leading] * pulsechargelist[following] > 0.0:
                    pulsechargelist[leading] += pulsechargelist[following]
                    pulsechargelist[following] = 0.0
                elif pulsechargelist[following] > 0.0:
                    leading = following
                following += 1
        else:
            # needs to be better
            for i in range(1, len(pulsetimelist)):
                if (
                    pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                ) < mingap and pulsechargelist[i] * pulsechargelist[i - 1] > 0.0:
                    mingap = pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                    minindex = i
            # If less than limit, combine pulses
            while mingap <= self.minTsep:
                if pulsechargelist[minindex] > pulsechargelist[minindex - 1]:
                    pulsechargelist[minindex] += pulsechargelist[minindex - 1]
                    pulsechargelist[minindex - 1] = 0.0
                else:
                    pulsechargelist[minindex - 1] += pulsechargelist[minindex]
                    pulsechargelist[minindex] = 0.0
                mingap = self.minTsep + 1.0
                minindex = -1
                # reestablish new min gap
                for i in range(1, len(pulsetimelist)):
                    if (
                        pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                    ) < mingap and pulsechargelist[i] * pulsechargelist[i - 1] > 0.0:
                        mingap = pulsetimelist[i][0] - pulsetimelist[i - 1][0]
                        minindex = i

        pmtinlist = []
        for i in range(len(pulsetimelist)):
            if pulsechargelist[-1 - i] < self.PEthreshold:
                continue
            if pulsetimelist[i][1] not in pmtinlist:
                pmtinlist.append(pulsetimelist[i][1])
        pmtinlist.sort()

        for i in range(len(pmtinlist)):
            newomkey = AddPMTKey(omkey, int(pmtinlist[i]))
            outputpulsemap[newomkey] = dataclasses.I3RecoPulseSeries()

        for i in range(len(pulsetimelist)):
            # remove pulses with too low charge.
            if pulsechargelist[-1 - i] < self.PEthreshold:
                continue
            rpulse = dataclasses.I3RecoPulse()
            rpulse.time = pulsetimelist[i][0]
            # saturate pulses with too much charge.
            if pulsechargelist[-1 - i] > self.PEsaturation:
                rpulse.charge = self.PEsaturation
                rpulse.charge += (pulsechargelist[-1 - i] - self.PEsaturation) * (
                    self.PEsaturation / pulsechargelist[-1 - i]
                )
            else:
                rpulse.charge = pulsechargelist[-1 - i]
            rpulse.width = int(pulsetimelist[i][1])
            if not (DOMpulsemap is None):
                if omkey not in DOMpulsemap.keys():
                    DOMpulsemap[omkey] = dataclasses.I3RecoPulseSeries()
                DOMpulsemap[omkey].append(rpulse)
            newomkey = AddPMTKey(omkey, pulsetimelist[i][1])
            outputpulsemap[newomkey].append(rpulse)

    # Apply the responce of the PMT, this includes pulse combining for pulses too close together.
    def ApplyPMTResponse(self, mcpemap, darkhits):
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        outputpulsemap_nonoise = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsemap = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsetimelist = {}

        self.npmt = 0

        for istring in range(1, self.nstring + 1):
            for iom in range(1, self.nom + 1):
                omkey = OMKey(istring, iom, 0)
                if omkey in mcpemap.keys():
                    RealPulses = self.ApplyPMTCharacteristics(mcpemap[omkey])
                    # print("real pulse = "+str(len(mcpemap[omkey])))
                    pulsetimelist = []
                    realpulsetimelist = []
                    isreal = []
                    DarkPulses = []
                    if omkey in darkhits.keys():
                        DarkPulses = self.ApplyPMTCharacteristicsDarkNoise(
                            darkhits[omkey]
                        )
                        # print("dark")
                        # print(len(DarkPulses))
                        # print(len(darkhits[omkey]))
                        j = 0
                        i = 0
                        while j < len(RealPulses) and i < len(DarkPulses):
                            if RealPulses[j][0] <= DarkPulses[i][0]:
                                pulsetimelist.append(RealPulses[j])
                                isreal.append(1.0)
                                j += 1
                            else:
                                pulsetimelist.append(DarkPulses[i])
                                isreal.append(0.0)
                                i += 1
                        while j < len(RealPulses):
                            pulsetimelist.append(RealPulses[j])
                            isreal.append(1.0)
                            j += 1
                        while i < len(DarkPulses):
                            pulsetimelist.append(DarkPulses[i])
                            isreal.append(0.0)
                            i += 1
                    else:
                        pulsetimelist = RealPulses
                        for i in range(len(pulsetimelist)):
                            isreal.append(1.0)

                    pulsechargelist = []
                    realpulsechargelist = []

                    for i in range(len(pulsetimelist)):
                        charge = self.randomService.gaus(
                            self.chargemean, self.chargesigma
                        )
                        pulsechargelist.append(charge)
                        if isreal[i] > 0.0:
                            realpulsechargelist.append(charge)

                    # print(len(RealPulses))
                    # print(len(realpulsechargelist))
                    # print(len(pulsetimelist))
                    # print(len(pulsechargelist))
                    # print("MakeRecos")

                    if len(DarkPulses) > 0 and len(RealPulses) > 0:
                        self.MakeRecoPulses(
                            RealPulses,
                            realpulsechargelist,
                            omkey,
                            outputpulsemap_nonoise,
                        )
                        self.MakeRecoPulses(
                            pulsetimelist,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )
                    elif len(RealPulses) > 0:
                        self.MakeRecoPulses(
                            RealPulses,
                            realpulsechargelist,
                            omkey,
                            outputpulsemap_nonoise,
                            DOMpulsemap,
                        )
                        for ipmt in range(1, self.dom_properties.GetNPMTs() + 1):
                            newomkey = AddPMTKey(omkey, ipmt)
                            if newomkey in outputpulsemap_nonoise.keys():
                                outputpulsemap[newomkey] = outputpulsemap_nonoise[
                                    newomkey
                                ]
                    elif len(DarkPulses) > 0:
                        self.MakeRecoPulses(
                            pulsetimelist,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )
                    else:
                        continue

                elif omkey in darkhits.keys():
                    if omkey in mcpemap.keys():
                        continue
                    if len(darkhits[omkey]) > 0:
                        DarkPulses = self.ApplyPMTCharacteristicsDarkNoise(
                            darkhits[omkey]
                        )
                        pulsechargelist = []
                        for i in range(len(DarkPulses)):
                            pulsechargelist.append(
                                self.randomService.gaus(
                                    self.chargemean, self.chargesigma
                                )
                            )

                        self.MakeRecoPulses(
                            DarkPulses,
                            pulsechargelist,
                            omkey,
                            outputpulsemap,
                            DOMpulsemap,
                        )

                # if omkey in darkhits.keys() or omkey in mcpemap.keys() :
                #    for ipmt in range(self.dom_properties.GetNPMTs()) :
                #        newomkey = AddPMTKey(omkey,ipmt)
                #        print(newomkey)
                #        print("all")
                #        if  newomkey in outputpulsemap.keys() :
                #            print(outputpulsemap[newomkey])
                #        else:
                #            print("no pulses")
                #        print("no noise")
                #        if newomkey in outputpulsemap_nonoise.keys() :
                #            print(outputpulsemap_nonoise[newomkey])
                #        else :
                #            print("no pulses")

        return outputpulsemap, outputpulsemap_nonoise, DOMpulsemap

    def Geometry(self, frame):
        if len(self.dropstrings) > 0:
            self.domsUsed = dataclasses.I3OMGeoMap()
            for omkey in frame["I3Geometry"].omgeo.keys():
                if omkey.string in self.dropstrings:
                    continue
                self.domsUsed[omkey] = frame["I3Geometry"].omgeo[omkey]
        else:
            self.domsUsed = frame["I3Geometry"].omgeo

        domkeylist = []
        self.nstring = 0
        self.nom = 0
        self.npmt = 0
        for omkey in self.domsUsed.keys():
            self.nstring = max(self.nstring, omkey.string)
            self.nom = max(self.nom, omkey.om)
            self.npmt = max(self.npmt, omkey.pmt)
            domkeylist.append(NoPMTKey(omkey))
        self.domkeys = set(domkeylist)
        self.domkeys = sorted(self.domkeys, key=lambda x: (x.string, x.om, x.pmt))
        self.nstring += 1
        self.nom += 1
        self.npmt += 1
        self.modgeo = frame["I3ModuleGeoMap"]
        self.PushFrame(frame)

    #########################################################################

    def Simulation(self, frame):
        if len(self.dropstrings) > 0:
            frame["SimulatedDOMs"] = self.domsUsed

        self.WroteActiveDOMsToSimFrame = True
        self.PushFrame(frame)

    def DAQ(self, frame):
        if not self.WroteActiveDOMsToSimFrame:
            simframe = icetray.I3Frame("S")
            self.Simulation(simframe)

        photonmap = frame[self.inputmap]
        noisepulses = [
            frame[noisepulsename]
            for noisepulsename in self.noisePulseSeries
            if frame.Has(noisepulsename)
        ]
        lennoisepulses = sum([len(noisepulses[i]) for i in range(len(noisepulses))])
        # print("ndoms = "+str(len(photonmap)))
        # split photons from DOMs to PMTs-on-DOMs
        # and apply weights
        if (len(photonmap) < 1) and self.noPureNoiseEvents:
            return
        if len(photonmap) < 1 and lennoisepulses < 1:
            return

        photonmap_on_pmts, frame[self.inputmap + "_pmtsplit"] = self.SplitPMTs(
            photonmap, self.dropstrings
        )
        # print(" pmt split pulses = " + str(len(photonmap_on_pmts)))
        # find the minimum and maximum time (with padding)
        max_pt, min_pt = self.GetMaxMinTimes(photonmap_on_pmts)
        # print("get min max")
        darkhits = {}
        if self.add_noise:
            # add dark noise
            # print("add noise")
            darkhits = self.AddDarkHits(max_pt, min_pt)
            darkhits = self.AddEnvironmentNoise(darkhits, noisepulses)
        # print("ndarkhitpmts = " + str(len(darkhits)))
        # Save the MCPEs to the frame
        # frame[self.outputmap_mcpe] = mcpemap
        # print("apply pmt response")
        (
            frame[self.outputmap],
            frame[self.outputmap + "_nonoise"],
            frame["triggerpulsemap"],
        ) = self.ApplyPMTResponse(photonmap_on_pmts, darkhits)
        # print(frame['triggerpulsemap'])
        # print("final")
        # print(frame[self.outputmap])

        self.PushFrame(frame)



class K40DOMSimulation(icetray.I3ConditionalModule):
    '''
    Simple implementation of the PMT response.

    This version uses a slightly more detailed POM
    acceptance model initially used for K40 studies
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('inputmap',
                          'Name of the I3Photons from clsim',
                          'I3Photons')
        self.AddParameter('outputmap',
                          'Name of the output I3RecoPulseSeriesMap',
                          'I3RecoPulseSeriesMap')
        self.AddParameter('outputmap_mcpe',
                          'Name of the output I3MCPESeriesMap',
                          'I3MCPESeriesMap')
        self.AddParameter('add_noise',
                          'Should random noise be added?',
                          True)
        self.AddParameter('PMT_tts',
                          'Transit time spread of PMT',
                          3.0 * I3Units.ns)
        self.AddParameter('PMT_ts',
                          'Transit time of PMT',
                          25.0 * I3Units.ns)
        self.AddParameter('chargesigma',
                          'Sigma of charge distribution',
                          0.3)
        self.AddParameter('chargemean',
                          'Mean of Charge distribution',
                          1.0)
        self.AddParameter('DNprob',
                          'Dark Noise rate (pulses per ns)',
                          0.000001)
        self.AddParameter('APprob',
                          'Total AP probability',
                          0.06)
        self.AddParameter('APmeantime_1',
                          'Mean of early time AP distribution',
                          2000.0 * I3Units.ns)
        self.AddParameter('APtimesigma_1',
                          'Sigma of early time AP distribution',
                          1000.0 * I3Units.ns)
        self.AddParameter('APmeantime_2',
                          'Mean of late time AP distribution',
                          8000.0 * I3Units.ns)
        self.AddParameter('APtimesigma_2',
                          'Sigma of late time AP distribution',
                          2000.0 * I3Units.ns)
        self.AddParameter('APComponetRatio',
                          'Fraction of AP in early time component',
                          0.3)
        self.AddParameter('LPprob',
                          'Probability for a late pulse',
                          0.01)
        self.AddParameter('minTsep',
                          ' Minimum time for a separated pulse.',
                          3.0)
        self.AddParameter('PEthreshold',
                          ' Pulse charge threshold',
                          0.25)
        self.AddParameter('PEsaturation',
                          'Saturation threshold for PMT',
                          100.0)
        self.AddParameter('RandomService',
                          'Random Service')
        self.AddParameter('SplitDoms',
                          '',
                          True)
        self.AddParameter('DOMAcceptanceFile',
                          '',
                          '')
        self.AddParameter('PMTQEFile',
                          '',
                          '')
        self.AddParameter('DropStrings',
                          '',
                          [])
        self.AddParameter('NoisePulseSeries',
                          '',
                          [])
        self.AddParameter('NoPureNoiseEvents',
                          '',
                          True)
        self.AddParameter('drop_empty',
                          'Bool to determine if empty frames should be removed from the i3 file',
                          False)
        self.AppParameter('use_dark',
                          'Bool to determine if dark noise will be included',
                          False)
        self.AppParameter('use_k40',
                          'Bool to determine if  noik40se will be included',
                          False)
        self.AppParameter('k40_map',
                          'Name of the k40 hit map',
                          'K40Hits')
        self.AddOutBox('OutBox')


    def Configure(self):
        self.inputmap          = self.GetParameter('inputmap')
        self.outputmap         = self.GetParameter('outputmap')
        self.outputmap_mcpe    = self.GetParameter('outputmap_mcpe')
        self.add_noise         = self.GetParameter('add_noise')
        self.PMT_tts           = self.GetParameter('PMT_tts')
        self.PMT_ts            = self.GetParameter('PMT_ts')
        self.chargesigma       = self.GetParameter('chargesigma')
        self.chargemean        = self.GetParameter('chargemean')
        self.DNprob            = self.GetParameter('DNprob')
        self.APprob            = self.GetParameter('APprob')
        self.APmeantime_1      = self.GetParameter('APmeantime_1')
        self.APtimesigma_1     = self.GetParameter('APtimesigma_1')
        self.APmeantime_2      = self.GetParameter('APmeantime_2')
        self.APtimesigma_2     = self.GetParameter('APtimesigma_2')
        self.APComponetRatio   = self.GetParameter('APComponetRatio')
        self.LPprob            = self.GetParameter('LPprob')
        self.minTsep           = self.GetParameter('minTsep')
        self.PEthreshold       = self.GetParameter('PEthreshold')
        self.PEsaturation      = self.GetParameter('PEsaturation')
        self.randomService     = self.GetParameter('RandomService')
        self.splitDOMs         = self.GetParameter('SplitDoms')
        self.dropstrings       = self.GetParameter('DropStrings')
        self.noisePulseSeries  = self.GetParameter('NoisePulseSeries')
        self.noPureNoiseEvents = self.GetParameter('NoPureNoiseEvents')
        self.drop_empty        = self.GetParameter('drop_empty')
        self.use_dark          = self.GetParameter('use_dark')
        self.use_k40           = self.GetParameter('use_k40')
        self.dark_map          = self.GetParameter('dark_map')
        self.k40_map           = self.GetParameter('k40_map')

        if self.GetParameter('PMTQEFile') != '':
            GetPMTQETable(self.photonweights, self.GetParameter('PMTQEFile'))

        kwargs_DOMProperties = {}
        if self.GetParameter('DOMAcceptanceFile') != '':
            kwargs_DOMProperties['PMTAcceptanceFile'] = self.GetParameter(
                'DOMAcceptanceFile'
            )
        if self.GetParameter('PMTQEFile') != '':
            kwargs_DOMProperties['PMTQEFile'] = self.GetParameter('PMTQEFile')

        # load the DOM properties from their respective configuration files
        self.dom_properties            = DOMProperties(**kwargs_DOMProperties)
        self.WroteActiveDOMsToSimFrame = False


        # load the updated module acceptance
        self.module = POM()
    

    def get_mcpe_map(self, pulse_map, drop_strings=[]):
        '''
        Read a split pmt pulse map from the frame and
        return an OM wide mcpe map of hit times and PMTs
        '''
        mcpe_map = {}
        
        # make new map with individual PMTs
        for pmtkey in pulse_map.keys():
            # ignore this omkey if it is supposed to be dropped
            if pmtkey.string in drop_strings:
                continue

            omkey = NoPMTKey(ModuleKey(pmtkey.string, pmtkey.om))
            if omkey not in mcpe_map.keys():
                mcpe_map[omkey] = []

            for pulse in pulse_map[pmtkey]:
                mcpe_map[omkey].append((pulse.time, pmtkey.pmt)) # mcpe map entries are tuples (time, pmt)

        return mcpe_map


    def add_environment_noise(self, dark_hits, noise_pulses):
        '''
        Adds environmental noise ??????????????????
        '''
        new_mcpe_map = {}
        merged_map   = {}

        for pulse_series in noise_pulses:
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeries()):
                for dom in pulse_series.keys():
                    nopmtkey = NoPMTKey(dom)
                    if nopmtkey not in new_mcpe_map.keys():
                        new_mcpe_map[nopmtkey] = list()
                    new_mcpe_map[nopmtkey].extend(
                        [
                            (pulse_series[dom][i].time, dom.pmt)
                            for i in range(len(pulse_series[dom]))
                        ]
                    )
            elif type(pulse_series) == type(simclasses, dataclasses.I3PhotonSeries()):
                for dom in noise_pulses.keys():
                    new_mcpe_map[nopmtkey(dom)] = list()
                    for pulse in noise_pulses[dom]:
                        pmtid = self.GetPMT(
                            photonDir=[pulse.dir.x, pulse.dir.y, pulse.dir.z],
                            wl=pulse.wavelength / I3Units.nanometer,
                            weight=pulse.weight,
                        )
                        new_mcpe_map[nopmtkey].append((pulse.time, pmtid))
            else:
                print('Invalid noise pulse series')
        for dom in new_mcpe_map.keys():
            new_mcpe_map[dom].sort(key=lambda x: x[0])

        for dom in new_mcpe_map.keys():
            merged_map[dom] = list()
            if dom in dark_hits.keys():
                i = 0
                j = 0
                while i < len(dark_hits[dom]) and j < len(new_mcpe_map[dom]):
                    if dark_hits[dom].time < sortedpulses[j].time:
                        merged_map[dom].append(dark_hits[dom])
                        i += 1
                    else:
                        merged_map[dom].append(sortedpulses[j])
                        j += 1
                if i < len(dark_hits[dom]):
                    merged_map[dom].append(dark_hits[dom])
                    i += 1
                if j < len(sortedpulses):
                    merged_map[dom].append(sortedpulses[j])
                    j += 1
            else:
                merged_map[dom] = new_mcpe_map[dom]

        for dom in dark_hits.keys():
            if dom not in merged_map.keys():
                merged_map[dom] = dark_hits[dom]
        return merged_map


    def combine_ordered_lists(self, list_1, order_1, list_2, order_2):
        '''
        Combines two ALREADY ORDERED lists into a single
        ordered list preserving the total order. Pass the
        lists with the objects you want to order along with
        arrays of what they need to be ordered by.

        Also returns an array of the same size as the final
        combined list filled with 1s and 2s corresponding
        to which list each element in the final list came from
        '''
        list_1_len = len(list_1)
        list_2_len = len(list_2)

        combined_list = np.empty(list_1_len + list_2_len, dtype=object)
        source_list   = np.zeros(list_1_len + list_2_len)

        i = 0
        j = 0

        # combine the first entries of the lists
        # based on comparing their times to eachother
        # to determine the order
        while i < list_1_len and j < list_2_len:
            if order_1[i] <= order_2[j]:
                combined_list[i + j] = list_1[i]
                source_list[i + j]   = 1
                i += 1
            else:
                combined_list[i + j] = list_2[j]
                source_list[i + j]   = 2
                j += 1

        # add the rest of list 1 if needed
        while i < list_1_len:
            combined_list[i + j] = list_1[i]
            source_list[i + j]   = 1
            i += 1

        # add the rest of list 2 if needed
        while j < list_2_len:
            combined_list[i + j] = list_2[j]
            source_list[i + j]   = 2
            j += 1

        return combined_list, source_list


    def apply_pmt_timing_characteristics(self, mcpe_map, late_pulses=True, after_pulses=True):
        '''
        Applies transit time and transit time spread
        to a given mcpe map

        late_pulses and after_pulses are two booleans
        which toggle including late pulses and afterpulses
        '''
        pulse_time_list = []
        # pulseseries   = dataclasses.I3RecoPulseSeries()
        # dompulseseries = dataclasses.I3RecoPulseSeries()
        pulse_out_of_order_time_list = []

        # for every photoelectron shift the time
        # based on tts and late pulse probability
        for pe in mcpe_map:
            # time shifted by tts
            time = self.randomService.gaus(pe[0], self.PMT_tts)

            # check if the pulse will be a late pulse based
            # on the late pulse probability
            if late_pulses:
                if self.randomService.uniform(0.0, 1.0) < self.LPprob:
                    time += self.randomService.gaus(self.PMT_ts * 2.0, np.sqrt(2.0) * self.PMT_tts)

            if len(pulse_time_list) < 1 or time > pulse_time_list[-1][0]:
                pulse_time_list.append((time, pe[1]))
            else:
                pulse_out_of_order_time_list.append((time, pe[1]))

        # if there are pulses out of order we need to combine the
        # two ordered lists into one
        if len(pulse_out_of_order_time_list) > 0:
            pulse_out_of_order_time_list.sort(key=lambda x: x[0])

            pulse_times              = [pe[0] for pe in pulse_time_list]
            pulse_out_of_order_times = [pe[0] for pe in pulse_out_of_order_time_list]

            pulse_time_list, _ = self.combine_ordered_lists(list_1  = pulse_time_list.copy(),
                                                            order_1 = pulse_times,
                                                            list_2  = pulse_out_of_order_time_list,
                                                            order_2 = pulse_out_of_order_times)

        # now add afterpulses
        if after_pulses:
            for pe in pulse_time_list:
                pulse_out_of_order_time_list = []
                if self.randomService.uniform(0.0, 1.0) < self.APprob:
                    if self.randomService.uniform(0.0, 1.0) < self.APComponetRatio:
                        time = pe[0] + self.randomService.gaus(self.APmeantime_1, self.APtimesigma_1)
                    else:
                        time = pe[0] + self.randomService.gaus(self.APmeantime_2, self.APtimesigma_2)
                    pulse_out_of_order_time_list.append((time, pe[1]))

            # if there are pulses or afteruplses out of order we need
            # to combine the two ordered lists into one
            if len(pulse_out_of_order_time_list) > 0:
                pulse_out_of_order_time_list.sort(key=lambda x: x[0])

                pulse_times              = [pe[0] for pe in pulse_time_list]
                pulse_out_of_order_times = [pe[0] for pe in pulse_out_of_order_time_list]

                pulse_time_list, _ = self.combine_ordered_lists(list_1  = pulse_time_list.copy(),
                                                                order_1 = pulse_times,
                                                                list_2  = pulse_out_of_order_time_list,
                                                                order_2 = pulse_out_of_order_times)

        return pulse_time_list


    def make_reco_pulse(self, pulse_time_list, pulse_charge_list, omkey, output_pulse_map, om_pulse_map=None):
        '''
        Populates the output_pulse_map with an I3RecoPulseSeries
        based on the input pulse times and charges
        '''
        min_gap   = 4.0
        min_index = -1

        if len(pulse_time_list) > 100:
            leading   = 0
            following = 1
            while following < len(pulse_time_list):
                if (
                    pulse_time_list[following][0] - pulse_time_list[leading][0]
                ) < 3.0 and pulse_charge_list[leading] * pulse_charge_list[following] > 0.0:
                    pulse_charge_list[leading] += pulse_charge_list[following]
                    pulse_charge_list[following] = 0.0
                elif pulse_charge_list[following] > 0.0:
                    leading = following
                following += 1
        else:
            # needs to be better
            for i in range(1, len(pulse_time_list)):
                if (
                    pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                ) < min_gap and pulse_charge_list[i] * pulse_charge_list[i - 1] > 0.0:
                    min_gap = pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                    min_index = i
            # If less than limit, combine pulses
            while min_gap <= self.minTsep:
                if pulse_charge_list[min_index] > pulse_charge_list[min_index - 1]:
                    pulse_charge_list[min_index] += pulse_charge_list[min_index - 1]
                    pulse_charge_list[min_index - 1] = 0.0
                else:
                    pulse_charge_list[min_index - 1] += pulse_charge_list[min_index]
                    pulse_charge_list[min_index] = 0.0
                min_gap = self.minTsep + 1.0
                min_index = -1
                # reestablish new min gap
                for i in range(1, len(pulse_time_list)):
                    if (
                        pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                    ) < min_gap and pulse_charge_list[i] * pulse_charge_list[i - 1] > 0.0:
                        min_gap = pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                        min_index = i

        pmt_in_list = []
        for i in range(len(pulse_time_list)):
            if pulse_charge_list[-1 - i] < self.PEthreshold:
                continue
            if pulse_time_list[i][1] not in pmt_in_list:
                pmt_in_list.append(pulse_time_list[i][1])
        pmt_in_list.sort()

        for i in range(len(pmt_in_list)):
            newomkey = AddPMTKey(omkey, pmt_in_list[i])
            output_pulse_map[newomkey] = dataclasses.I3RecoPulseSeries()

        for i in range(len(pulse_time_list)):
            # remove pulses with too low charge.
            if pulse_charge_list[-1 - i] < self.PEthreshold:
                continue
            rpulse = dataclasses.I3RecoPulse()
            rpulse.time = pulse_time_list[i][0]
            # saturate pulses with too much charge.
            if pulse_charge_list[-1 - i] > self.PEsaturation:
                rpulse.charge = self.PEsaturation
                rpulse.charge += (pulse_charge_list[-1 - i] - self.PEsaturation) * (
                    self.PEsaturation / pulse_charge_list[-1 - i]
                )
            else:
                rpulse.charge = pulse_charge_list[-1 - i]
            rpulse.width = pulse_time_list[i][1]
            if not (om_pulse_map is None):
                if omkey not in om_pulse_map.keys():
                    om_pulse_map[omkey] = dataclasses.I3RecoPulseSeries()
                om_pulse_map[omkey].append(rpulse)
            newomkey = AddPMTKey(omkey, pulse_time_list[i][1])
            output_pulse_map[newomkey].append(rpulse)


    def OLDapply_pmt_response(self, mcpe_map, noise_mcpe_map):
        '''
        Apply the response of the PMT, including combining pulses
        that are too close together
        '''
        output_pulse_map         = dataclasses.I3RecoPulseSeriesMap()
        output_pulse_map_nonoise = dataclasses.I3RecoPulseSeriesMap()
        om_pulse_map             = dataclasses.I3RecoPulseSeriesMap()
        # DOMpulsetimelist = {}

        # self.npmt = 0

        for string_index in range(1, self.nstring + 1):
            for om_index in range(1, self.nom + 1):
                omkey = OMKey(string_index, om_index, 0)

                # if there was a real hit at this omkey
                if omkey in mcpe_map.keys():
                    real_pulses = self.apply_pmt_timing_characteristics(mcpe_map[omkey])

                    pulse_time_list      = []
                    real_pulse_time_list = []

                    is_real     = []
                    noise_pulses = []
                    if omkey in noise_mcpe_map.keys():
                        noise_pulses = self.apply_pmt_timing_characteristics(noise_mcpe_map[omkey], late_pulses=False)

                        # combine the two pulse time lists keeping
                        # track of which pulse is real and which
                        # is noise
                        pulse_times = [pe[0] for pe in real_pulses]
                        noise_times  = [pe[0] for pe in noise_pulses]
                        
                        pulse_time_list, source = self.combine_ordered_lists(list_1  = real_pulses,
                                                                             order_1 = pulse_times,
                                                                             list_2  = noise_pulses,
                                                                             order_2 = noise_times)
                        
                        is_real = np.abs(source - 2)

                    # if there are no dark hits just mark
                    # every hit as being real
                    else:
                        pulse_time_list = real_pulses
                        is_real         = np.ones(len(pulse_time_list))

                    pulse_charge_list      = []
                    real_pulse_charge_list = []

                    # collect the charges associated with
                    # each pulse time and them to the real
                    # charge list if they are not from noise
                    for i in range(len(pulse_time_list)):
                        charge = self.randomService.gaus(self.chargemean, self.chargesigma)
                        pulse_charge_list.append(charge)
                        if is_real[i] > 0.0:
                            real_pulse_charge_list.append(charge)

                    # turn lists of pulse times and charges into reco pulses

                    # if there are both dark pulses and real pulses
                    if len(noise_pulses) > 0 and len(real_pulses) > 0:
                        self.make_reco_pulse(real_pulses, real_pulse_charge_list, omkey,
                                             output_pulse_map_nonoise)
                        self.make_reco_pulse(pulse_time_list, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)
                    
                    # if there are only real pulses
                    elif len(real_pulses) > 0:
                        self.make_reco_pulse(real_pulses, real_pulse_charge_list, omkey,
                                             output_pulse_map_nonoise, om_pulse_map)
                        
                        # add this to the regular output noise map as well
                        for pmt_index in range(0, self.dom_properties.GetNPMTs()):
                            new_omkey = AddPMTKey(omkey, pmt_index)
                            if new_omkey in output_pulse_map_nonoise.keys():
                                output_pulse_map[new_omkey] = output_pulse_map_nonoise[new_omkey]
                    
                    # if there are only noise pulses
                    elif len(noise_pulses) > 0:
                        self.make_reco_pulse(pulse_time_list, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)
                    else:
                        continue

                # if there was only a dark hit at this om
                elif omkey in noise_mcpe_map.keys():
                    # avoid double counting
                    if omkey in mcpe_map.keys():
                        continue

                    if len(noise_mcpe_map[omkey]) > 0:
                        noise_pulses = self.apply_pmt_timing_characteristics(noise_mcpe_map[omkey], late_pulses=False)

                        # for every dark pulse add a charge
                        pulse_charge_list = []
                        for i in range(len(noise_pulses)):
                            pulse_charge_list.append(self.randomService.gaus(self.chargemean, self.chargesigma))

                        # make reco pulses and add to the output map
                        self.make_reco_pulse(noise_pulses, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)

        return output_pulse_map, output_pulse_map_nonoise, om_pulse_map


    def apply_pmt_response(self, mcpe_map):
        '''
        Apply the response of the PMT, including combining pulses
        that are too close together
        '''
        mcpe_map = self.apply_dead_time(mcpe_map.copy)

        output_pulse_map = dataclasses.I3RecoPulseSeriesMap()
        om_pulse_map     = dataclasses.I3RecoPulseSeriesMap()

        for omkey in mcpe_map.keys():
            pulse_charge_list = []

            # collect the charges associated with
            # each pulse time and them to the real
            # charge list if they are not from noise
            for i in range(len(mcpe_map)):
                charge = self.randomService.gaus(self.chargemean, self.chargesigma)
                pulse_charge_list.append(charge)
            
            self.make_reco_pulse(mcpe_map, pulse_charge_list, omkey, output_pulse_map, om_pulse_map)
        
        return output_pulse_map, om_pulse_map



        output_pulse_map_nonoise = dataclasses.I3RecoPulseSeriesMap()
        om_pulse_map             = dataclasses.I3RecoPulseSeriesMap()
        # DOMpulsetimelist = {}

        # self.npmt = 0

        for string_index in range(1, self.nstring + 1):
            for om_index in range(1, self.nom + 1):
                omkey = OMKey(string_index, om_index, 0)

                # if there was a real hit at this omkey
                if omkey in mcpe_map.keys():
                    real_pulses = self.apply_pmt_timing_characteristics(mcpe_map[omkey])

                    pulse_time_list      = []
                    real_pulse_time_list = []

                    is_real     = []
                    noise_pulses = []
                    if omkey in noise_mcpe_map.keys():
                        noise_pulses = self.apply_pmt_timing_characteristics(noise_mcpe_map[omkey], late_pulses=False)

                        # combine the two pulse time lists keeping
                        # track of which pulse is real and which
                        # is noise
                        pulse_times = [pe[0] for pe in real_pulses]
                        noise_times  = [pe[0] for pe in noise_pulses]
                        
                        pulse_time_list, source = self.combine_ordered_lists(list_1  = real_pulses,
                                                                             order_1 = pulse_times,
                                                                             list_2  = noise_pulses,
                                                                             order_2 = noise_times)
                        
                        is_real = np.abs(source - 2)

                    # if there are no dark hits just mark
                    # every hit as being real
                    else:
                        pulse_time_list = real_pulses
                        is_real         = np.ones(len(pulse_time_list))

                    pulse_charge_list      = []
                    real_pulse_charge_list = []

                    # collect the charges associated with
                    # each pulse time and them to the real
                    # charge list if they are not from noise
                    for i in range(len(pulse_time_list)):
                        charge = self.randomService.gaus(self.chargemean, self.chargesigma)
                        pulse_charge_list.append(charge)
                        if is_real[i] > 0.0:
                            real_pulse_charge_list.append(charge)

                    # turn lists of pulse times and charges into reco pulses

                    # if there are both dark pulses and real pulses
                    if len(noise_pulses) > 0 and len(real_pulses) > 0:
                        self.make_reco_pulse(real_pulses, real_pulse_charge_list, omkey,
                                             output_pulse_map_nonoise)
                        self.make_reco_pulse(pulse_time_list, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)
                    
                    # if there are only real pulses
                    elif len(real_pulses) > 0:
                        self.make_reco_pulse(real_pulses, real_pulse_charge_list, omkey,
                                             output_pulse_map_nonoise, om_pulse_map)
                        
                        # add this to the regular output noise map as well
                        for pmt_index in range(0, self.dom_properties.GetNPMTs()):
                            new_omkey = AddPMTKey(omkey, pmt_index)
                            if new_omkey in output_pulse_map_nonoise.keys():
                                output_pulse_map[new_omkey] = output_pulse_map_nonoise[new_omkey]
                    
                    # if there are only noise pulses
                    elif len(noise_pulses) > 0:
                        self.make_reco_pulse(pulse_time_list, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)
                    else:
                        continue

                # if there was only a dark hit at this om
                elif omkey in noise_mcpe_map.keys():
                    # avoid double counting
                    if omkey in mcpe_map.keys():
                        continue

                    if len(noise_mcpe_map[omkey]) > 0:
                        noise_pulses = self.apply_pmt_timing_characteristics(noise_mcpe_map[omkey], late_pulses=False)

                        # for every dark pulse add a charge
                        pulse_charge_list = []
                        for i in range(len(noise_pulses)):
                            pulse_charge_list.append(self.randomService.gaus(self.chargemean, self.chargesigma))

                        # make reco pulses and add to the output map
                        self.make_reco_pulse(noise_pulses, pulse_charge_list, omkey,
                                             output_pulse_map, om_pulse_map)

        return output_pulse_map, output_pulse_map_nonoise, om_pulse_map


    def apply_dead_time(self, mcpe_map):
        '''
        Removes successive hits on the same PMT
        based on the PMT dead time

        !!! ignores dark hits for now but these should
        !!! also be included
        '''
        dead_time             = 10.
        dead_removed_mcpe_map = {}

        for omkey in mcpe_map.keys():
            dead_removed_mcpe_map[omkey] = []

            last_hit_time = -9999.
            for pe in mcpe_map[omkey]:
                if pe[0] - last_hit_time > dead_time:
                    dead_removed_mcpe_map[omkey].append((pe[0], pe[1]))
                    last_hit_time = pe[0]
        
        return dead_removed_mcpe_map


    def merge_mcpe_maps(self, mcpe_map_1, mcpe_map_2):
        '''
        Combines two mcpe maps into one
        keeping the mcpes time ordered
        '''
        merged_map = {}

        # find the unique omkeys
        all_omkeys = np.unique(list(mcpe_map_1.keys()) + list(mcpe_map_2.keys()))

        for omkey in all_omkeys:
            # if the omkey is in one of the maps no need to go through the whole merge
            if omkey not in mcpe_map_1.keys():
                merged_map[omkey] = mcpe_map_2[omkey]
            if omkey not in mcpe_map_2.keys():
                merged_map[omkey] = mcpe_map_1[omkey]
            
            times_1 = [pe[0] for pe in mcpe_map_1[omkey]]
            times_2 = [pe[0] for pe in mcpe_map_2[omkey]]

            combined_mcpes, _ = self.combine_ordered_lists(list_1  = mcpe_map_1[omkey],
                                                           order_1 = times_1,
                                                           list_2  = mcpe_map_2[omkey],
                                                           order_2 = times_2)
            merged_map[omkey] = combined_mcpes
        
        return merged_map


    def Geometry(self, frame):
        # filter out the strings we want to use based
        # on dropstrings
        if len(self.dropstrings) > 0:
            self.domsUsed = dataclasses.I3OMGeoMap()
            for omkey in frame['I3Geometry'].omgeo.keys():
                if omkey.string in self.dropstrings:
                    continue
                self.domsUsed[omkey] = frame['I3Geometry'].omgeo[omkey]
        else:
            self.domsUsed = frame['I3Geometry'].omgeo


        domkeylist   = []
        self.nstring = 0
        self.nom     = 0
        self.npmt    = 0

        for omkey in self.domsUsed.keys():
            self.nstring = max(self.nstring, omkey.string)
            self.nom     = max(self.nom, omkey.om)
            self.npmt    = max(self.npmt, omkey.pmt)
            domkeylist.append(NoPMTKey(omkey))

        self.domkeys  = set(domkeylist)
        self.domkeys  = sorted(self.domkeys, key=lambda x: (x.string, x.om, x.pmt))
        self.nstring += 1
        self.nom     += 1
        self.npmt    += 1

        self.PushFrame(frame)


    def Simulation(self, frame):
        if len(self.dropstrings) > 0:
            frame['SimulatedDOMs'] = self.domsUsed

        self.WroteActiveDOMsToSimFrame = True

        self.PushFrame(frame)


    def DAQ(self, frame):
        if not self.WroteActiveDOMsToSimFrame:
            simframe = icetray.I3Frame('S')
            self.Simulation(simframe)
        
        simulation_pulse_map = frame[self.inputmap]
        simulation_mcpe_map  = self.get_mcpe_map(simulation_pulse_map, self.dropstrings)

        length_noise_pulses = 0
        if self.use_dark:
            dark_pulse_map       = frame[self.dark_map]
            dark_mcpe_map        = self.get_mcpe_map(dark_pulse_map, self.dropstrings)
            length_noise_pulses += len(dark_pulse_map)
        if self.use_k40:
            k40_pulse_map        = frame[self.k40_map]
            k40_mcpe_map         = self.get_mcpe_map(k40_pulse_map, self.dropstrings)
            length_noise_pulses += len(k40_pulse_map)
        
        # if there are no pulses or noise
        # just drop the frame
        if self.drop_empty:
            if (len(simulation_pulse_map) < 1) and self.noPureNoiseEvents:
                return            
            if len(simulation_pulse_map) < 1 and length_noise_pulses < 1:
                return

        noise_mcpe_maps = []

        #all_omkeys = []
        # apply the pmt timing characteristics to each mcpe map
        for omkey in simulation_mcpe_map.keys():
            simulation_mcpe_map[omkey] = self.apply_pmt_timing_characteristics(simulation_mcpe_map[omkey].copy)
            # if omkey not in all_omkeys:
            #     all_omkeys.append(omkey)
        
        if self.use_dark:
            for omkey in dark_mcpe_map.keys():
                dark_mcpe_map[omkey] = self.apply_pmt_timing_characteristics(dark_mcpe_map[omkey].copy, late_pulses=False)
                # if omkey not in all_omkeys:
                #     all_omkeys.append(omkey)
            noise_mcpe_maps.append(dark_mcpe_map)
        
        if self.use_k40:
            for omkey in k40_mcpe_map.keys():
                k40_mcpe_map[omkey] = self.apply_pmt_timing_characteristics(k40_mcpe_map[omkey].copy)
                # if omkey not in all_omkeys:
                #     all_omkeys.append(omkey)
            noise_mcpe_maps.append(k40_mcpe_map)

        noise_mcpe_map = {}
        # merge all the noise
        for i in len(noise_mcpe_maps):
            noise_mcpe_map = self.merge_mcpe_maps(noise_mcpe_map, noise_mcpe_maps[i])

        # merge the signal and noise
        total_mcpe_map = self.merge_mcpe_maps(simulation_mcpe_map, noise_mcpe_map)

        # apply the PMT pulse response without any noise
        (no_noise_output_pulses, _) = self.apply_pmt_response(simulation_mcpe_map)

        # apply the PMT pulse response to the combined signal and noise
        (output_pulses, om_pulses) = self.apply_pmt_response(total_mcpe_map)
        




        # combine all noise mcpe maps into one
        
        # dark_hits = {}
        # if self.add_noise:
        #     max_pt, min_pt = self.get_dark_noise_time_bounds(simulation_mcpe_map)
        #     dark_hits      = self.generate_dark_hits(max_pt, min_pt)
        #     #dark_hits      = self.add_environment_noise(dark_hits, noise_pulses) # not defined right now (commented out)

        # # account for PMT dead time DO THIS AFTER TTS!!!
        # dead_removed_photon_map = self.apply_dead_time(simulation_mcpe_map)

        # (output_pulses, output_pulses_nonoise, om_pulses,) = self.apply_pmt_response(dead_removed_photon_map, dark_hits)

        frame[self.outputmap]              = output_pulses
        frame[self.outputmap + '_nonoise'] = no_noise_output_pulses
        frame['triggerpulsemap']           = om_pulses

        self.PushFrame(frame)