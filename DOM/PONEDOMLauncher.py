from icecube import icetray
from icecube import dataclasses
from icecube.icetray import I3Units

class SimpleDOMSimulation(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile","GCD to be simulated",'')
        self.AddParameter("inputmap","Name of the Physics I3MCTree name","I3MCPulseSeriesMap")
        self.AddParameter("outputmap","Name of the noise I3MCTree name","I3RecoPulseSeriesMap")
        self.AddParameter("PMT_tts","Noise Hits for mDOMs to be injected",1.0)
        self.AddParameter("PMT_ts","Noise Hits for mDOMs to be injected",25.0)
        self.AddParameter("chargesigma","Noise Hits for mDOMs to be injected",0.5)
        self.AddParameter("chargemean","Noise Hits for mDOMs to be injected",1.0)
        self.AddParameter("DNprob"," for mDOMs to be injected",0.000001)
        self.AddParameter("APprob"," for mDOMs to be injected",0.15)
        self.AddParameter("APmeantime_1"," for mDOMs to be injected",2000.)
        self.AddParameter("APtimesigma_1"," for mDOMs to be injected",1000.)
        self.AddParameter("APmeantime_2"," for mDOMs to be injected",8000.)
        self.AddParameter("APtimesigma_2"," for mDOMs to be injected",2000.)
        self.AddParameter("APComponetRatio"," for mDOMs to be injected",0.3)
        self.AddParameter("LPprob"," for mDOMs to be injected",0.01)
        self.AddParameter("minTsep"," for mDOMs to be injected",10.0)
        self.AddParameter("PEthreshold"," for mDOMs to be injected",0.25)
        self.AddParameter("PEsaturation"," for mDOMs to be injected",100.0)

        self.AddOutBox("OutBox")

    def Configure(self):

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

	def DAQ(self,frame) :

		outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
		mcpulsemap = frame[self.inputmap]
		pulsetimelist = []
		sqrt2 = 1.414213562373095

		for omkey, pulses in mcpulsemap:
			recopulsemap[omkey] = dataclasses.I3RecoPulseSeries()
			#trueHists
			i=0
			for pulse in pulses:
				time = random_service.gaus(pulse.time,self.PMT_tts*I3Units.ns)
				if random_service.uniform() < self.LPprob :
					time += random_service.gaus(self.PMT_ts,sqrt2*self.PMT_tts*I3Units.ns)
				pulsetimelist.append(time)
				
			#darkhits
			max_pt = max(pulsetimelist) - 1000.
			min_pt = min(pulsetimelist) + 1000.
			ndarkhists = random_service.poisson((max_pt-min_pt)*self.DNprob)
			for i in range(ndarkhists) :
				time = random_service.uniform(min_pt,max_pt)
				pulsetimelist.append(time)

			for time in pulsetimelist :
				if random_service.uniform() < self.APprob :
					if random_service.uniform() < self.APComponetRatio :
						time = time + random_service.gaus(self.APmeantime_1*I3Units.ns,self.APtimesigma_1*I3Units.ns)
						pulsetimelist.append(time)
					else :
						time = time + random_service.gaus(self.APmeantime_2*I3Units.ns,self.APtimesigma_2*I3Units.ns)
						pulsetimelist.append(time)
			
			if len(pulsetimelist)<1 :
				frame[self.outputmap] = outputpulsemap
				return

			#time  order pulse
			sort(pulsetimelist)
			pulsechargelist = []
			#combine pulses that are too close
			for i in range(pulsetimelist):
				pulsechargelist.append(random_service.gaus(chargemean,chargesigma))
				j = i-1
				k = i+1
			while k<(len(pulsetimelist)-1) and pulsetimelist[k]-pulsetimelist[k-1] < minTsep :
				pulsechargelist[i] += random_service.gaus(chargemean,chargesigma)
				pulsechargelist.append(0.0)
				k+=1
			i=k-1

			for i in range(pulsetimelist):
				#remove pulses with too low charge.
				if pulsechargelist[i]<self.PEthreshold :
					continue
				rpulse = dataclasses.I3RecoPulse()
				rpulse.time = pulsetimelist[i]
				#saturate pulses with too much charge.
				if pulsechargelist[i] > self.PEsaturation :
					rpulse.charge = self.PEsaturation
					rpulse.charge += (pulsechargelist[i]-self.PEsaturation)*(self.PEsaturation/pulsechargelist[i])
				else :
					rpulse.charge = pulsechargelist[i]
				outputpulsemap[omkey].append(rpulse)

		frame[self.outputmap] = outputpulsemap