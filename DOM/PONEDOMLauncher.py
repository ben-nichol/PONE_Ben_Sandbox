from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import GetPMTAcceptance, GetPMTQETable, GetPMTQE, GetMaxTotalAcceptance, GetMaxAngularAcceptance

#split MCPhoton hits into PMTs on the DOMs.
def SplitPMTs(mcpulsemap,random_service) :
    newmcpulsemap = {}
    #make new map with individual PMTs
    total=0
    passed = 0
    for omkey in mcpulsemap.keys():
        for pulse in mcpulsemap[omkey]:
            pmtid = GetPMT(pulse.dir,pulse.wavelength,random_service.uniform(0.0,1.0))
            if pmtid < 0 :
                continue
            mcpulse = simclasses.I3Photon()
            mcpulse.time = pulse.time
            mcpulse.dir = pulse.dir
            mcpulse.pos = pulse.pos

            newomkey = OMKey(omkey.string, omkey.om, pmtid)
            if newomkey in newmcpulsemap.keys() :
                newmcpulsemap[newomkey].append(mcpulse)
            else :
                newmcpulsemap[newomkey] = simclasses.I3PhotonSeries()
                newmcpulsemap[newomkey].append(mcpulse)
    return  newmcpulsemap

#Get the min and max pulse times to set time window for dark hits.
def GetMaxMinTimes(mcpulsemap) :
    max_pt = -999999999.
    min_pt = 999999999.
    for omkey in mcpulsemap.keys():
        for pulse in mcpulsemap[omkey]:
            if max_pt < pulse.time :
                max_pt = pulse.time
            if min_pt > pulse.time :
                min_pt = pulse.time

    max_pt += 1000.
    min_pt -= 1000.

    return max_pt, min_pt

#Add Dark hits across all DOMs
def AddDarkHits(domsUsed,mcpulsemap,random_service,DNprob,max_pt,min_pt,splitDOMs) :
    for omkey in domsUsed.keys() :
        npmts = 1
        if splitDOMs :
            npmts = GetNPMTs()
        for ipmt in range(npmts) :
            key = OMKey(omkey.string,omkey.om,ipmt)
            ndarkhists = random_service.poisson((max_pt-min_pt)*DNprob)
            if ndarkhists < 1 :
                continue
            if key not in mcpulseOMkeys :
                mcpulsemap[key] = simclasses.I3PhotonSeries()
            for i in range(ndarkhists) :
                time = random_service.uniform(min_pt,max_pt)
                mcpulsemap[key].append(time*I3Units.ns)

#Apply the responce of the PMT, this includes pulse combining for pulses too close together. 
def ApplyPMTResponce(mcpulsemap,random_service,PMT_tts,PMT_ts,LPprob,APprob,APComponetRatio,APmeantime_1,APmeantime_2,APtimesigma_1,APtimesigma_2,chargemean,chargesigma,PEsaturation,PEthreshold) :
    outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
    outputmcpulsemap = simclasses.I3MCPulseSeriesMap()


    for omkey in mcpulsemap.keys():
        pulsetimelist = []
        pulseseries = dataclasses.I3RecoPulseSeries()
        mcpulseseries = simclasses.I3MCPulseSeries()
        #trueHists
        i=0
        for pulse in mcpulsemap[omkey]:
            time = random_service.gaus(pulse.time,PMT_tts*I3Units.ns)
            if random_service.uniform(0.0,1.0) < LPprob :
                time += random_service.gaus(self.PMT_ts,np.sqrt(2.0)*PMT_tts*I3Units.ns)
            pulsetimelist.append(time)

        for time in pulsetimelist :
            if random_service.uniform(0.0,1.0) < APprob :
                if random_service.uniform(0.0,1.0) < APComponetRatio :
                    time = time + random_service.gaus(APmeantime_1*I3Units.ns,APtimesigma_1*I3Units.ns)
                else :
                    time +=  random_service.gaus(APmeantime_2*I3Units.ns,APtimesigma_2*I3Units.ns)
                pulsetimelist.append(time)

        pulsetimelist.sort()
        pulsechargelist = []

        for i in range(len(pulsetimelist)) :
            mcpulse = simclasses.I3MCPulse()
            mcpulse.time = pulsetimelist[i]
            mcpulseseries.append(mcpulse)
            pulsechargelist.append(1.0)

        mingap = 4.0
        minindex = -1
        #Get MinGap
        for i in range(1,len(pulsetimelist)) :
            if (pulsetimelist[i]-pulsetimelist[i-1]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
                mingap = (pulsetimelist[i]-pulsetimelist[i-1])
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
                if (pulsetimelist[i]-pulsetimelist[i-1]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
                    mingap = (pulsetimelist[i]-pulsetimelist[i-1])
                    minindex = i

        for i in range(len(pulsechargelist)) :
            if pulsechargelist[i]>0.0 :
                pulsechargelist[i] = random_service.gaus(chargemean*pulsechargelist[i],np.sqrt(pulsechargelist[i])*chargesigma)

        for i in range(len(pulsetimelist)):
            #remove pulses with too low charge.
            if pulsechargelist[-1-i]<PEthreshold :
                continue
            rpulse = dataclasses.I3RecoPulse()
            rpulse.time = pulsetimelist[i]
            #saturate pulses with too much charge.
            if pulsechargelist[-1-i] > PEsaturation :
                rpulse.charge = PEsaturation
                rpulse.charge += (pulsechargelist[-1-i]-PEsaturation)*(PEsaturation/pulsechargelist[-1-i])
            else :
                rpulse.charge = pulsechargelist[-1-i]
            pulseseries.append(rpulse)
        newomkey = OMKey(omkey.string, omkey.om, omkey.pmt)

        outputpulsemap[newomkey]=pulseseries
        outputmcpulsemap[newomkey] = mcpulseseries

    return outputpulsemap, outputmcpulsemap


class SimpleDOMSimulation(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("inputmap","Name of the Physics I3MCTree name","I3MCPulseSeriesMap")
        self.AddParameter("outputmap","Name of the noise I3MCTree name","I3RecoPulseSeriesMap")
        self.AddParameter("PMT_tts","Transit time spread of PMT",1.0)
        self.AddParameter("PMT_ts","Transit time of PMT",25.0)
        self.AddParameter("chargesigma","Sigma of charge distribution",0.3)
        self.AddParameter("chargemean","Mean of Charge distribution ",1.0)
        self.AddParameter("DNprob","Dark Noise rate (pulses per ns)",0.000001)
        self.AddParameter("APprob","Total AP probability",0.06)
        self.AddParameter("APmeantime_1","Mean of early time AP distribution",2000.)
        self.AddParameter("APtimesigma_1","Sigma of early time AP distribution",1000.)
        self.AddParameter("APmeantime_2"," mean of late time AP distribution",8000.)
        self.AddParameter("APtimesigma_2"," sigma of late time AP distribution",2000.)
        self.AddParameter("APComponetRatio","fraction of AP in early time component",0.3)
        self.AddParameter("LPprob","Probability for a late pulse",0.01)
        self.AddParameter("minTsep"," minimum time for a separated pulse.",3.0)
        self.AddParameter("PEthreshold"," Pulse charge threshold",0.25)
        self.AddParameter("PEsaturation"," Saturation threshold for PMT",100.0)
        self.AddParameter("RandomService","Random Service")
        self.AddParameter("GenWaveforms","Generate Waveforms?",False)
        self.AddParameter("SplitDoms","",True)
        self.AddParameter("DOMAcceptanceFile","","")
        self.AddParameter("PMTQEFile","","")
        self.AddParameter("QEBaseValue","",0.4)
        self.AddParameter("AcceptBaseValue","",1.0)
        self.AddOutBox("OutBox")

    def Configure(self):
        global QEBase
        global AccBase

        self.gcdFile = self.GetParameter("GCDFile")
        self.inputmap = self.GetParameter("inputmap")
        self.outputmap = self.GetParameter("outputmap")
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
        self.genWaveforms = self.GetParameter("GenWaveforms")
        self.splitDOMs = self.GetParameter("SplitDoms")
        QEBase = self.GetParameter("QEBaseValue")
        if self.GetParameter("DOMAcceptanceFile") != "" :
            GetPMTAcceptanceTable(self.GetParameter("DOMAcceptanceFile"))
            if self.GetParameter("AcceptBaseValue") < 0.0 :
                AccBase = GetMaxAngularAcceptance()
            else :
                AccBase = self.GetParameter("AcceptBaseValue")
        else :
            self.splitDOMs = False
        if self.GetParameter("PMTQEFile") != "" :
            GetPMTQETable(self.GetParameter("PMTQEFile"))

    def DAQ(self,frame) :

        random_service = self.randomService
        outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
        outputmcpulsemap = simclasses.I3MCPulseSeriesMap()
        mcpulsemap = frame[self.inputmap]
        mcpulseOMKeys = mcpulsemap.keys()

        domsUsed = frame['I3Geometry'].omgeo

        # print("splitting pulses")
        if self.splitDOMs :
            mcpulsemap = SplitPMTs(mcpulsemap,random_service)
            mcpulseOMKeys = mcpulsemap.keys()

        max_pt, min_pt = GetMaxMinTimes(mcpulsemap)

        #add dark noise

        AddDarkHits(domsUsed,mcpulsemap,random_service,self.DNprob,max_pt,min_pt,self.splitDOMs)

        frame[self.outputmap], frame[self.outputmap+"_MCpulses"] = ApplyPMTResponce(mcpulsemap,
                                                            random_service,
                                                            self.PMT_tts,
                                                            self.PMT_ts,
                                                            self.LPprob,
                                                            self.APprob,
                                                            self.APComponetRatio,
                                                            self.APmeantime_1,
                                                            self.APmeantime_2,
                                                            self.APtimesigma_1,
                                                            self.APtimesigma_2,
                                                            self.chargemean,
                                                            self.chargesigma,
                                                            self.PEsaturation,
                                                            self.PEthreshold)

        self.PushFrame(frame)

