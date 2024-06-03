from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, DOMProperties, Geant4PMTAcceptance

# from DOMUtility import NoPMTKey, AddPMTKey, DOMProperties


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

            pmt_ids = self.pmt_acc.check_pmt_hit(photon_positions, photon_wavelengths, photon_weights)

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