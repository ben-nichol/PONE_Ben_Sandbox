from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from DOM.pulseshape import specharge, pulsewave



def GenerateWaveform(time,pulsetimes,pulsecharge,mv_to_adc):

  V = np.random.normal(0.0,1.0)
  for i in range(len(pulsetimes)) :
    V += pulsecharge[i]*pulsewave(time-pulsetimes[i])*mv_to_adc

  return int(V)

class WaveformBuilder(icetray.I3ConditionalModule):
    """
    Simple Implementation of the PMT response.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("inputmap","Name of the Physics I3MCTree name","I3MCPulseSeriesMap")
        self.AddParameter("outputmap","Name of the noise I3MCTree name","I3RecoPulseSeriesMap")
    
        self.AddOutBox("OutBox")

    def Configure(self):
        self.inputmap = self.GetParameter("inputmap")   
        self.outputmap = self.GetParameter("outputmap")
        self.mv_to_ADC = 10.
        self.SPECharge = self.mv_to_ADC

    def DAQ(self,frame) :
    
        outputWaveformmap = dataclasses.I3WaveformSeriesMap()   
        mcpulsemap = frame[self.inputmap]
        mcpulseOMKeys = mcpulsemap.keys()
        sqrt2 = 1.414213562373095

        for omkey in mcpulseOMKeys:
            waveformseries = dataclasses.I3WaveformSeries()
     
            waveform = []
            pulsetimes = []
            pulsecharge = []
            charge = []
      
            for pulse in mcpulsemap[omkey]:
                pulsetimes.append(pulse.time)
                pulsecharge.append(pulse.charge)
            binwidth = 4.0
            tmin = min(pulsetimes)-300.
            tmax = max(pulsetimes)+300.
            nbins = int((tmax-tmin)/binwidth)

            #assume 3 bin noise
            for i in range(nbins) :
                waveform.append(GenerateWaveform(float(int(tmin+i*binwidth)),pulsetimes,pulsecharge,self.mv_to_ADC))

            splitwaveform = []
            binstart = []
            binstop = []
            inwaveform = False

            for i in range(len(waveform)-10) :
                if inwaveform :
                    if abs(waveform[i]) < 4 and abs(waveform[i+1]) < 4 and abs(waveform[i+2]) < 4 and abs(waveform[i+3]) < 4:
                        for j in range(10) :
                            splitwaveform[-1].append(waveform[i+j])
                        binstop.append(i+j)
                        inwaveform = False
                    else:
                        splitwaveform[-1].append(waveform[i])
                if abs(waveform[i]) > 4 :
                    splitwaveform.append([])
                    inwaveform = True
                    binstart.append(max(i-9,0))
                    for j in range(binstart[-1],i+1) :
                        splitwaveform[-1].append(waveform[j])

            for i in range(len(splitwaveform)) :
                i3waveform = dataclasses.I3Waveform()
                i3waveform.waveform = splitwaveform[i]
                i3waveform.time = binstart[i]*binwidth+tmin
                i3waveform.binwidth = binwidth
                waveformseries.append(i3waveform)

            #print(waveformseries)
            outputWaveformmap[omkey] = waveformseries


#
#      deriv = [0.0]
#      for i in range(1,len(waveform)) :
#        deriv.append(waveform[i]-waveform[i-1]);
#
#      #pulsefinding
#      Derthreshold = 6
#      StopDerivthresh = 3
#      StopVoltthresh = 3
#      baseDerthresh = 3
#      baseVolthresh = 5
#      DefaultBaseline = 2.0
#      PulseSplitVoltThresh = -6
#      baselinesamples = []
#      pulsefindertimes = []
#      pulsefindercharge = []
#      pulsefindermaxderiv = []
#      pulsefinderstartbin = []
#      pulsefinderstopbin = []
#      pulsefindermaxV = []
#      pulsefindermavbin = []
#      inpulse = False
#      cansplit = False
#
#      localbaseline = DefaultBaseline
#      for i in range(len(waveform)) :
#        if inpulse :
#
#          if deriv[i] < PulseSplitVoltThresh :
#            cansplit = True
#
#          if deriv[i] < StopDerivthresh and waveform[i] < StopVoltthresh :
#            pulsefinderstopbin.append(i)
#            inpulse = False
#            cansplit = False
#
#          if cansplit and deriv[i] > Derthreshold :
#            pulsefinderstopbin.append(i)
#            pulsefindercharge.append(0.0)
#            pulsefindermaxderiv.append(deriv[i])
#            pulsefinderstartbin.append(i)
#            cansplit = False
#
#          if deriv[i]>pulsefindermaxderiv[-1]:
#            pulsefindermaxderiv[-1] = deriv[i]
#
#          pulsefindercharge[-1] += float(waveform[i])-localbaseline
#          
#        elif deriv[i] > Derthreshold :
#          inpulse = True
#          if len(baselinesamples) < 5 :
#            localbaseline = sum(baselinesamples)/len(baselinesamples)
#          pulsefindercharge.append(float(waveform[i])-localbaseline)
#          pulsefindermaxderiv.append(deriv[i])
#          pulsefinderstartbin.append(i)
#
#        elif waveform[i] < baseVolthresh and deriv[i] < baseDerthresh :
#          baselinesamples.append(waveform[i])

#      newomkey = OMKey(omkey.string, omkey.om, 0)
#      for i in range(len(pulsefindertimes)):
#        rpulse = dataclasses.I3RecoPulse()
#        rpulse.time = pulsefindertimes[i]
#        rpulse.charge = pulsefindercharge[i]/self.SPECharge
#        rpulse.width = float(pulsefinderstopbin[i]-pulsefinderstartbin[i])*3.0
#        pulseseries.append(rpulse)
#
#      outputpulsemap[newomkey]= pulseseries 
#
#
        frame[self.outputmap] = outputWaveformmap
    
        self.PushFrame(frame) 
