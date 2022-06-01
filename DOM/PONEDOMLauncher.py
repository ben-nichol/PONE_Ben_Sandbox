from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, DOMProperties

class SimpleDOMSimulation(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("inputmap","Name of the I3Photons from clsim","I3Photons")
        self.AddParameter("outputmap","Name of the output I3RecoPulseSeriesMap","I3RecoPulseSeriesMap")
        self.AddParameter("outputmap_mcpe","Name of the output I3MCPESeriesMap","I3MCPESeriesMap")
        self.AddParameter("add_noise","Should random noise be added?",True)
        self.AddParameter("PMT_tts","Transit time spread of PMT",1.0*I3Units.ns)
        self.AddParameter("PMT_ts","Transit time of PMT",25.0*I3Units.ns)
        self.AddParameter("chargesigma","Sigma of charge distribution",0.3)
        self.AddParameter("chargemean","Mean of Charge distribution ",1.0)
        self.AddParameter("DNprob","Dark Noise rate (pulses per ns)",0.000001)
        self.AddParameter("APprob","Total AP probability",0.06)
        self.AddParameter("APmeantime_1","Mean of early time AP distribution",2000.*I3Units.ns)
        self.AddParameter("APtimesigma_1","Sigma of early time AP distribution",1000.*I3Units.ns)
        self.AddParameter("APmeantime_2"," mean of late time AP distribution",8000.*I3Units.ns)
        self.AddParameter("APtimesigma_2"," sigma of late time AP distribution",2000.*I3Units.ns)
        self.AddParameter("APComponetRatio","fraction of AP in early time component",0.3)
        self.AddParameter("LPprob","Probability for a late pulse",0.01)
        self.AddParameter("minTsep"," minimum time for a separated pulse.",3.0)
        self.AddParameter("PEthreshold"," Pulse charge threshold",0.25)
        self.AddParameter("PEsaturation"," Saturation threshold for PMT",100.0)
        self.AddParameter("RandomService","Random Service")
        self.AddParameter("SplitDoms","",True)
        self.AddParameter("DOMAcceptanceFile","","")
        self.AddParameter("PMTQEFile","","")
        self.AddParameter("DropStrings","",[])
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
        
        if self.GetParameter("PMTQEFile") != "" :
            GetPMTQETable(self.photonweights,self.GetParameter("PMTQEFile"))

        kwargs_DOMProperties = {}
        if self.GetParameter("DOMAcceptanceFile") != "" :
            kwargs_DOMProperties["PMTAcceptanceFile"] = self.GetParameter("DOMAcceptanceFile")
        if self.GetParameter("PMTQEFile") != "" :
            kwargs_DOMProperties["PMTQEFile"] = self.GetParameter("PMTQEFile")

        # load the DOM properties from their respective configuration files
        self.dom_properties = DOMProperties(**kwargs_DOMProperties)

    #########################################################################

    #split Photon hits into PMTs on the DOMs.
    def SplitPMTs(self, photonmap, dropstrings=[]) :
        newphotonmap = {}
        #make new map with individual PMTs

        for omkey in photonmap.keys():
            if omkey.string in dropstrings :
                continue
            newomkey = NoPMTKey(omkey)
            newphotonmap[newomkey] = []
            for pulse in photonmap[omkey]:
                pmtid = self.GetPMT(
                    photonDir = [pulse.dir.x,pulse.dir.y,pulse.dir.z],
                    wl = pulse.wavelength/I3Units.nanometer,
                    weight = pulse.weight
                    )

                if pmtid < 0 :
                    continue
                newphotonmap[newomkey].append((pulse.time,pmtid))
        #print("split")
        #print(newphotonmap)
        return newphotonmap

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
    def GetPMT(self, photonDir, wl,weight):
        random = self.randomService.uniform(0.0,1.0)

        theta = 0.0
        phi = 0.0
        if len(photonDir) < 3 :
            theta = photonDir[0]
            phi = photonDir[1]
        else :
            theta = np.arccos(photonDir[2])
            phi = np.arctan2(photonDir[1],photonDir[0])

        while phi < 0.0 :
            phi += 2.*np.pi

        thetaBin = max(0,min(178,int(180.0*theta/np.pi)))
        phiBin = max(0,min(358,int(180.0*phi/np.pi)))
        pmtprobs = []
        for i in range(len(self.dom_properties.PMTacceptance)) :
            pmtprobs.append(self.dom_properties.PMTacceptance[i][thetaBin][phiBin]/self.dom_properties.maxAngularAcceptance)

        totalprob = sum(pmtprobs)
        if random > totalprob :
            #print("angular kill photon")
            return -1

        random = self.randomService.uniform(0.0,1.0)
        i=0
        sumprob = pmtprobs[0]
        while random > sumprob/totalprob and i<len(pmtprobs)-1:
            i+=1
            sumprob += pmtprobs[i]

        qe = self.dom_properties.GetPMTQEnm(wl)
        # note: the probaility should only be weight * qe. It also includes maxAngularAcceptance
        #       here, since we already took this factor into account in `GetPMT` (which also rejects
        #       hits.)
        probability = weight * qe * self.dom_properties.maxAngularAcceptance

        #print("probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format( probability, weight, qe, self.dom_properties.maxAngularAcceptance, wl ))

        if probability > 1.:
            print("ERROR: probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format( probability, weight, qe, self.dom_properties.maxAngularAcceptance, wl ))
            raise RuntimeError("The combined detection probability should never be > 1. You need to re-generate your I3Photons with a higher overall bias weight (setting `WavelengthAcceptance`)")

        ### reject photon according to `probability`
        random = self.randomService.uniform(0.0,1.0)
        if random > probability:
            #print("kill photon")
            return -1
        
        #print("pmt = "+str(i+1))

        return i+1

    #Get the min and max pulse times to set time window for dark hits.
    def GetMaxMinTimes(self, mcpemap) :
        max_pt = -999999999.
        min_pt = 999999999.
        #print(mcpemap)
        for omkey in mcpemap.keys():
            for pulse in mcpemap[omkey]:
                if max_pt < pulse[0] :
                    max_pt = pulse[0]
                if min_pt > pulse[0] :
                    min_pt = pulse[0]

        max_pt += 10000.
        min_pt -= 1000.

        return max_pt, min_pt

    # Add Dark hits across all DOMs (addes I3MCPEs)
    def AddDarkHits(self,mcpemap, max_pt, min_pt):
        newmcpemap = {}
        #print("min = "+str(min_pt)+" max = "+str(max_pt))
        npmts = int(self.dom_properties.GetNPMTs())
        for omkey in self.domkeys :
            #print(omkey)
            dt = self.randomService.exp(1./(self.DNprob*npmts)) + min_pt
            ipmt = 1 + self.randomService.integer(npmts)
            if omkey not in mcpemap.keys() :
                #print("not in")
                first = True
                while dt < max_pt :
                    if first :
                        newmcpemap[omkey] = []
                        first = False
                    #print(dt)
                    newmcpemap[omkey].append((dt,ipmt))
                    dt += self.randomService.exp(1./(self.DNprob*npmts))
                    ipmt = 1 + self.randomService.integer(npmts)
            else :
                #print("in")
                newmcpemap[omkey] = []
                for pulse in mcpemap[omkey] :
                    while dt < pulse[0] :
                        newmcpemap[omkey].append((dt,ipmt))
                        dt += self.randomService.exp(1.0/(self.DNprob*npmts))
                        ipmt = 1 + self.randomService.integer(npmts)
                        #print(dt)
                    #print(pulse)
                    newmcpemap[omkey].append(pulse)
                while dt < max_pt :
                    newmcpemap[omkey].append((dt,ipmt))
                    dt += self.randomService.exp(1.0/(self.DNprob*npmts))
                    ipmt = 1 + self.randomService.integer(npmts)
                    #print(dt)
        #print("dark")
        #print(newmcpemap)

        return newmcpemap

    #Apply the responce of the PMT, this includes pulse combining for pulses too close together. 
    def ApplyPMTResponse(self, mcpemap) :
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsemap = dataclasses.I3RecoPulseSeriesMap()
        DOMpulsetimelist = {}
        for omkey in mcpemap.keys():
            pulsetimelist = []
            pulseseries = dataclasses.I3RecoPulseSeries()
            dompulseseries = dataclasses.I3RecoPulseSeries()
            pulseoutoftimelist = []
            for pe in mcpemap[omkey]:
                time = self.randomService.gaus(pe[0],self.PMT_tts)
                if self.randomService.uniform(0.0,1.0) < self.LPprob :
                    time += self.randomService.gaus(self.PMT_ts,np.sqrt(2.0)*self.PMT_tts)
                if len(pulsetimelist)<1 or time > pulsetimelist[-1][0] :
                    pulsetimelist.append((time,pe[1]))
                else :
                    pulseoutoftimelist.append((time,pe[1]))
            if len(pulseoutoftimelist) > 0 :
                intimelist = pulsetimelist.copy()
                pulsetimelist = []
                pulseoutoftimelist.sort(key=lambda x:x[0])
                j=0
                i=0
                while j<len(intimelist) and i<len(pulseoutoftimelist):
                    if intimelist[j][0] <= pulseoutoftimelist[i][0] :
                        pulsetimelist.append(intimelist[j])
                        j+=1
                    else :
                        pulsetimelist.append(pulseoutoftimelist[i])
                        i+=1

            for pe in pulsetimelist :
                pulseoutoftimelist = []
                if self.randomService.uniform(0.0,1.0) < self.APprob :
                    if self.randomService.uniform(0.0,1.0) < self.APComponetRatio :
                        time = pe[0] + self.randomService.gaus(self.APmeantime_1,self.APtimesigma_1)
                    else :
                        time = pe[0] + self.randomService.gaus(self.APmeantime_2,self.APtimesigma_2)
                    pulseoutoftimelist.append((time,pe[1]))

            if len(pulseoutoftimelist) > 0 :
                intimelist = pulsetimelist.copy()
                pulsetimelist = []
                pulseoutoftimelist.sort(key=lambda x:x[0])
                j=0
                i=0
                while j<len(intimelist) and i<len(pulseoutoftimelist):
                    if intimelist[j][0] <= pulseoutoftimelist[i][0] :
                        pulsetimelist.append(intimelist[j])
                        j+=1
                    else :
                        pulsetimelist.append(pulseoutoftimelist[i])
                        i+=1

            pulsechargelist = []

            for i in range(len(pulsetimelist)) :
                pulsechargelist.append(1.0)

            #splitPMTs

            mingap = 4.0
            minindex = -1
            #Get MinGap
            for i in range(1,len(pulsetimelist)) :
                if (pulsetimelist[i][0]-pulsetimelist[i-1][0]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
                    mingap = (pulsetimelist[i][0]-pulsetimelist[i-1][0])
                    minindex = i
            #If less than limit, combine pulses
            while mingap <= 3.0 :
                if pulsechargelist[minindex] > pulsechargelist[minindex-1]:
                    pulsechargelist[minindex] += pulsechargelist[minindex-1]
                    pulsechargelist[minindex-1] = 0.0
                else:
                    pulsechargelist[minindex-1] += pulsechargelist[minindex]
                    pulsechargelist[minindex] = 0.0
                mingap = 4.0
                minindex = -1
                #reestablish new min gap
                for i in range(1,len(pulsetimelist)) :
                    if (pulsetimelist[i][0]-pulsetimelist[i-1][0]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
                        mingap = (pulsetimelist[i][0]-pulsetimelist[i-1][0])
                        minindex = i

            for i in range(len(pulsechargelist)) :
                if pulsechargelist[i]>0.0 :
                    pulsechargelist[i] = self.randomService.gaus(self.chargemean*pulsechargelist[i],np.sqrt(pulsechargelist[i])*self.chargesigma)

            for i in range(len(pulsetimelist)):
                #remove pulses with too low charge.
                if pulsechargelist[-1-i]<self.PEthreshold :
                    continue
                rpulse = dataclasses.I3RecoPulse()
                rpulse.time = pulsetimelist[i][0]
                #saturate pulses with too much charge.
                if pulsechargelist[-1-i] > self.PEsaturation :
                    rpulse.charge = self.PEsaturation
                    rpulse.charge += (pulsechargelist[-1-i]-self.PEsaturation)*(self.PEsaturation/pulsechargelist[-1-i])
                else :
                    rpulse.charge = pulsechargelist[-1-i]
                rpulse.width = pulsetimelist[i][1]
                if omkey not in DOMpulsemap.keys() :
                    DOMpulsemap[omkey] = dataclasses.I3RecoPulseSeries()
                DOMpulsemap[omkey].append(rpulse)
                newomkey = AddPMTKey(omkey,pulsetimelist[i][1])
                if newomkey not in outputpulsemap.keys() :
                    outputpulsemap[newomkey] = dataclasses.I3RecoPulseSeries()
                outputpulsemap[newomkey].append(rpulse)
            
        return outputpulsemap, DOMpulsemap

    def Geometry(self,frame) :
        self.domsUsed = frame['I3Geometry'].omgeo
        domkeylist = []
        for omkey in self.domsUsed.keys() :
            domkeylist.append(NoPMTKey(omkey))
        self.domkeys = set(domkeylist)
        self.domkeys = sorted(self.domkeys, key = lambda x: (x.string, x.om))
        self.PushFrame(frame)
    #########################################################################
    def DAQ(self,frame) :
        #print("new event")
        photonmap = frame[self.inputmap]
        #print(photonmap)
        #print("get photons")
        # split photons from DOMs to PMTs-on-DOMs
        # and apply weights
        if len(photonmap) < 1 :
            return
        photonmap_on_pmts = self.SplitPMTs(photonmap, self.dropstrings)
        #print("split pulses")
        #frame["I3PhotonSplitPMT"] = mcpemap

        # find the minimum and maximum time (with padding)
        max_pt, min_pt = self.GetMaxMinTimes(photonmap_on_pmts)
        #print("get min max")
        if self.add_noise:
            # add dark noise
            #print("add noise")
            mcpemap = self.AddDarkHits(photonmap_on_pmts,max_pt,min_pt)

        # Save the MCPEs to the frame
        #frame[self.outputmap_mcpe] = mcpemap
        #print("apply pmt response")
        frame[self.outputmap], frame['triggerpulsemap'] = self.ApplyPMTResponse(mcpemap)

        #print("final")
        #print(frame[self.outputmap])

        self.PushFrame(frame)
