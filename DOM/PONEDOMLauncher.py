from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np

PMTacceptance = list()
PMTQE = list()
QEBase = 1.0
AccBase = 1.0

def GetPMTAcceptanceMax(infile) :

    domaccFile = open(infile,"r")
    lines = domaccFile.readlines()
    maxTotaleff = 0.0;

    PMTacceptance = list()

    zenithcount = 0
    for line in lines :
        splitline = line.split(" ",1000)
        if zenithcount % 179 == 0 :
            zenithcount = 0
            PMTacceptance.append([])
        PMTacceptance[-1].append([])
        for value in splitline :
            PMTacceptance[-1][-1].append(float(value))
        zenithcount += 1

    for theta in range(179) :
        for phi in range(359) :
            sum_acceptance = 0.0
            for ipmt in range(len(PMTacceptance)) :
                sum_acceptance += PMTacceptance[ipmt][theta][phi]
            if maxTotaleff < sum_acceptance :
                maxTotaleff = sum_acceptance

    return maxTotaleff

def GetPMTQETable(infile) :
    global PMTQE

    PMTQE_value = list()
    PMTQE_wl = list()

    pmtqeFile = open(infile,"r")
    lines = pmtqeFile.readlines()
    for line in lines :
        splitline = line.split(",",100)
        PMTQE_wl.append(float(splitline[0]))
        PMTQE_value.append(float(splitline[1]))
    #print("len(PMTQE_wl) = "+str(len(PMTQE_wl)))
    j = 1
    for i in range(1000):
     #   print("PMTQE_wl[0] = "+str(PMTQE_wl[0])+" PMTQE_wl[-1] = "+str(PMTQE_wl[-1]))
        if i<PMTQE_wl[0] :
            PMTQE.append(0.0)
            continue
        elif i>PMTQE_wl[-1] :
            PMTQE.append(0.0)
            continue
        
        while i > PMTQE_wl[j] :
            j += 1
        #print("j = "+str(j))
        #print("PMTQE_value[j-1] = "+str(PMTQE_value[j-1])+" PMTQE_value[j] = "+str(PMTQE_value[j]))
        #print("(float(i)-PMTQE_wl[j-1]) = "+str((float(i)-PMTQE_wl[j-1])))
        #print("(PMTQE_wl[j]-PMTQE_wl[j-1]) = "+str((PMTQE_wl[j]-PMTQE_wl[j-1])))
        QE = PMTQE_value[j-1] + (PMTQE_value[j]-PMTQE_value[j-1])*((float(i)-PMTQE_wl[j-1])/(PMTQE_wl[j]-PMTQE_wl[j-1]))
        PMTQE.append(QE)
        #print(" Wl = "+str(i)+" value ="+str(PMTQE[-1]))

def GetPMTQE(wl):
    global PMTQE


    #print("len(PMTQE)-1 = "+str(len(PMTQE)-1))

    if int(wl*1.0e9) > len(PMTQE)-1 :
        return 0.0

    return PMTQE[int(wl*1.0e9)]

def GetPMTAcceptanceTable(infile) :
    global PMTacceptance

    domaccFile = open(infile,"r")
    lines = domaccFile.readlines()

    zenithcount = 0
    for line in lines :

        splitline = line.split(" ",1000)
        if zenithcount % 179 == 0 :
            zenithcount = 0
            PMTacceptance.append([])
        PMTacceptance[-1].append([])
        for value in splitline :
            PMTacceptance[-1][-1].append(float(value))
        zenithcount += 1
            
def GetPMT(photonDir,wl,random):
    global PMTacceptance
    global QEBase
    global AccBase
    #print("wavelength = "+str(wl))
    #print("random = "+str(random))    
    QEProb = GetPMTQE(wl)/QEBase
    #print("QEProb = "+str(QEProb))
    if random > QEProb :
        return -1
  
    random = random/QEProb

    thetaBin = max(0,min(178,int(180.0*photonDir.zenith/np.pi)))
    phiBin = max(0,min(358,int(180.0*photonDir.azimuth/np.pi)))
    pmtprobs = []
    for i in range(len(PMTacceptance)) :
        pmtprobs.append(PMTacceptance[i][thetaBin][phiBin]/AccBase)

    #print("totalacc = "+str(sum(pmtprobs)))
    totalprob = sum(pmtprobs)
    if random > totalprob :
        return -1
    i=0;
    sumprob = pmtprobs[i]
    while random > sumprob/totalprob and i<len(pmtprobs)-1 :
       i+=1
       sumprob += pmtprobs[i]

    return i

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
    self.AddParameter("SplitDoms","",False)
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
            AccBase = GetPMTAcceptanceMax(self.GetParameter("DOMAcceptanceFile"))
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
    sqrt2 = 1.414213562373095

    if self.splitDOMs :
        newmcpulsemap = {}
        #make new map with individual PMTs
        total=0
        passed = 0
        for omkey in mcpulseOMKeys:
           for pulse in mcpulsemap[omkey]:
               pmtid = GetPMT(pulse.dir,pulse.wavelength,random_service.uniform(0.0,1.0))
               total+=1
               if pmtid < 0 :
                    continue
               passed +=1
               pmtid = 0
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
        mcpulsemap = newmcpulsemap
        mcpulseOMKeys = mcpulsemap.keys()
        #print(str(float(passed))+"/"+str(float(total)))

    for omkey in mcpulseOMKeys:
      pulsetimelist = []
      pulseseries = dataclasses.I3RecoPulseSeries()
      mcpulseseries = simclasses.I3MCPulseSeries()
      #trueHists
      i=0
      for pulse in mcpulsemap[omkey]:
        time = random_service.gaus(pulse.time,self.PMT_tts*I3Units.ns)
        if random_service.uniform(0.0,1.0) < self.LPprob :
          time += random_service.gaus(self.PMT_ts,sqrt2*self.PMT_tts*I3Units.ns)
        pulsetimelist.append(time)
      #darkhits
      max_pt = max(pulsetimelist) + 1000.
      min_pt = 0.0 #min(pulsetimelist) + 1000.
      
      ndarkhists = random_service.poisson((max_pt-min_pt)*self.DNprob)
      for i in range(ndarkhists) :
        time = random_service.uniform(min_pt,max_pt)
        pulsetimelist.append(time*I3Units.ns)
      nAP = 0
      for time in pulsetimelist :
        if random_service.uniform(0.0,1.0) < self.APprob :
          nAP+=1
          if random_service.uniform(0.0,1.0) < self.APComponetRatio :
            time = time + random_service.gaus(self.APmeantime_1*I3Units.ns,self.APtimesigma_1*I3Units.ns)
            pulsetimelist.append(time)
          else :
            time = time + random_service.gaus(self.APmeantime_2*I3Units.ns,self.APtimesigma_2*I3Units.ns)
            pulsetimelist.append(time)
      if len(pulsetimelist)<1 :
        frame[self.outputmap] = outputpulsemap
        return 
      #time  order pulse
      pulsetimelist.sort()
      pulsechargelist = []
      #combine pulses that are too close

      for i in range(len(pulsetimelist)) :
        mcpulse = simclasses.I3MCPulse()
        mcpulse.time = pulsetimelist[i]
        mcpulseseries.append(mcpulse)
        pulsechargelist.append(1.0)

      mingap = 4.0
      minindex = -1
      for i in range(1,len(pulsetimelist)) :
        if (pulsetimelist[i]-pulsetimelist[i-1]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
          mingap = (pulsetimelist[i]-pulsetimelist[i-1])
          minindex = i
      while mingap <= 3.0 : 
        if pulsechargelist[minindex] > pulsechargelist[minindex-1]:
          pulsechargelist[minindex] += pulsechargelist[minindex-1]
          pulsechargelist[minindex-1] = 0.0
        else:
          pulsechargelist[minindex-1] += pulsechargelist[minindex]
          pulsechargelist[minindex] = 0.0
        mingap = 4.0
        minindex = -1
        for i in range(1,len(pulsetimelist)) :
          if (pulsetimelist[i]-pulsetimelist[i-1]) < mingap and pulsechargelist[i]*pulsechargelist[i-1] > 0.0:
            mingap = (pulsetimelist[i]-pulsetimelist[i-1])
            minindex = i
        
      for i in range(len(pulsechargelist)) :
        if pulsechargelist[i]>0.0 :
          pulsechargelist[i] = random_service.gaus(self.chargemean*pulsechargelist[i],np.sqrt(pulsechargelist[i])*self.chargesigma)

      for i in range(len(pulsetimelist)):
        #remove pulses with too low charge.
        if pulsechargelist[-1-i]<self.PEthreshold :
          continue
        rpulse = dataclasses.I3RecoPulse()
        rpulse.time = pulsetimelist[i]
        #saturate pulses with too much charge.
        if pulsechargelist[-1-i] > self.PEsaturation :
          rpulse.charge = self.PEsaturation
          rpulse.charge += (pulsechargelist[-1-i]-self.PEsaturation)*(self.PEsaturation/pulsechargelist[-1-i])
        else :
          rpulse.charge = pulsechargelist[-1-i]
        pulseseries.append(rpulse)
      newomkey = OMKey(omkey.string, omkey.om, omkey.pmt)
      outputpulsemap[newomkey]=pulseseries 
      outputmcpulsemap[newomkey] = mcpulseseries  

    frame[self.outputmap] = outputpulsemap
    frame[self.outputmap+"_MCpulses"] = outputmcpulsemap
    
    self.PushFrame(frame) 

