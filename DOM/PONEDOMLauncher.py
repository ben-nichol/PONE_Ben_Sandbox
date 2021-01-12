from icecube import icetray
from icecube import dataclasses

def SimpleDOMSimulation(frame,inputmap="I3MCPulseSeriesMap",outputmap="I3RecoPulseSeriesMap",
						tts = 1.0,ts=25.0,chargesigma = 0.5,chargemean=1.0,DNprob=0.000001,
						APprob=0.15,APmeantime=3000.,APtimesigma=500.,LPprob=0.01,minTsep=10.0,
						PEthreshold=0.25,PEsaturation=100.0) :

	from icecube.icetray import I3Units

	outputpulsemap = dataclasses.I3RecoPulseSeriesMap()
	mcpulsemap = frame[inputmap]
	pulsetimelist = []
	sqrt2 = 1.414213562373095

	for omkey, pulses in mcpulsemap:
		recopulsemap[omkey] = dataclasses.I3RecoPulseSeries()
		#trueHists
		i=0
		for pulse in pulses:
			time = random_service.gaus(pulse.time,tts*I3Units.ns)
			if random_service.uniform() < LPprob :
				time += random_service.gaus(ts,sqrt2*tts*I3Units.ns)
			pulsetimelist.append(time)
			
			#add afterpulses
			newpulse = random_service.uniform() < APprob
			while newpulse :
				time = time + random_service.gaus(APmeantime*I3Units.ns,APtimesigma*I3Units.ns)
				pulsetimelist.append(time)
				newpulse = random_service.uniform() < APprob
				
		#darkhits
		max_pt = max(pulsetimelist) - 1000.
		min_pt = min(pulsetimelist) + 1000.
		ndarkhists = random_service.poisson((max_pt-min_pt)*DNprob)
		for i in range(ndarkhists) :
			time = random_service.uniform(min_pt,max_pt)
			pulsetimelist.append(time)

			newpulse = random_service.uniform() < APprob
			while newpulse :
				time = time + random_service.gaus(APmeantime*I3Units.ns,APtimesigma*I3Units.ns)
				pulsetimelist.append(time)
				newpulse = random_service.uniform() < APprob
			
		if len(pulsetimelist)<1 :
			frame[outputmap] = outputpulsemap
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
			if pulsechargelist[i]<PEthreshold :
				continue
			rpulse = dataclasses.I3RecoPulse()
			rpulse.time = pulsetimelist[i]
			#saturate pulses with too much charge.
			if pulsechargelist[i] > PEsaturation :
				rpulse.charge = PEsaturation
				rpulse.charge += (pulsechargelist[i]-PEsaturation)*(PEsaturation/pulsechargelist[i])
			else :
				rpulse.charge = pulsechargelist[i]
			outputpulsemap[omkey].append(rpulse)

	frame[outputmap] = outputpulsemap