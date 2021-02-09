#!/usr/bin/env python                                                                                                
 
# Import some useful ICECUBE modules                                                                                  
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
from scipy import special as sp                        # For the Gamma function 
import sys, os
from iminuit import Minuit
import argparse
import math as m

def LikelihoodFunctor(data,pdf,time_lim,dist_lim):
                
    pulse_series = data

    c = 0.299792458                                 # speed of light 
    n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/n                                       # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 557                                       # time parameter that has to be fit using simulations or data 
    pdf_tables = pdf
    table_time_lim = time_lim
    table_dist_lim = dist_lim
    darkprob = 1e-5

    def GetProbability(time,distance) :
  
      if time < table_time_lim[0] or time > table_time_lim[1] :
        return 0.0
      if distance < table_dist_lim[0] or distance > table_dist_lim[1] :
        return 0.0

      dist_bin = max(min(int(distance),len(pdf_tables)-1),0)
      time_bin = max(min(int(time-table_time_lim[0]),len(pdf_tables[dist_bin])-1),0)

      return pdf_tables[dist_bin][time_bin]

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(t0,d0,t1,d1,br): 
        darkrate = 1.e-7

        nloglike = 0.0

        for dom in pulse_series.keys() :

            for pulse in pulse_series[dom] :
                t = pulse.time
                charge = pulse.charge
                
                prob = br*GetProbability(t-t0,d0)
                prob += (1.-br)*GetProbability(t-t1,d1)
                prob += darkrate

                nloglike += -charge*np.log(prob)

        return nloglike

    return likelihoodFunction

class curveFit(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("pulseseries","Name of the Merged MCPE tree name","MergedSeriesMap")
        self.AddParameter("output","Track to store fit.","taufit")
        self.AddParameter("tables","tablesfile","")
        self.AddParameter("HitsInDOMsCut","Cut in the num fits in DOMS",200)
        self.AddOutBox("OutBox")

    def ReadTables(self) :

        infile = open(self.tablesfiles,"r")
        lines = infile.readlines()
        linecount = 0
        self.pdf = [[]]
        xcount = 0
        for line in lines :
            splitline = line.split(",",100)
            if linecount == 0 :
                nx = int(splitline[0].replace("\n",""))
                ny = int(splitline[1].replace("\n",""))
                minx = float(splitline[2].replace("\n",""))
                maxx = float(splitline[3].replace("\n",""))
                miny = float(splitline[4].replace("\n",""))
                maxy = float(splitline[5].replace("\n",""))
            else :
                if xcount == nx :
                    self.pdf.append([])
                    xcount = 0
                for value in splitline :
                    self.pdf[-1].append(float(value.replace("\n","")))
                    xcount += 1
            linecount += 1
        self.time_lim = [minx,maxx]
        self.dist_lim = [miny,maxy]

    def Configure(self):

        self.pulseseries = self.GetParameter("pulseseries")
        self.output = self.GetParameter("output")
        self.cuts = self.GetParameter("HitsInDOMsCut")

        self.c = 0.299792458                                 # speed of light 
        self.n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
        self.c_n = self.c/self.n                             # light in water
        self.theta_c = np.arccos(1./self.n)                  # Cherenkov angle in water in radians
        self.lambda_s = 120.                                 # scattering length of light for violet light
        self.lambda_a = 15.                                  # absorption length of light for violet light
        self.tau = 557                                       # time parameter that has to be fit using simulations or data   
 
        self.tablesfiles = self.GetParameter("Tables")
        if self.tablesfiles == "" :
          self.tablesfiles = default=os.getenv('PONESRCDIR')+"/data/fittertables.dat"

        self.ReadTables() 

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.                                   
    def DAQ(self,frame): 
        

        data = frame[self.pulseseries]

        biGauss_valuesMap = dataclasses.I3MapKeyVectorDouble()
        doublePeak_valuesMap = dataclasses.I3MapKeyVectorDouble()
        exitStatusMap = dataclasses.I3MapKeyVectorDouble()

        for omkey in data.keys():

            recoPulseList = data[omkey]
            recoPulse_timeList = np.array([recoPulse.time for recoPulse in recoPulseList])
            recoPulse_chargeList = np.array([recoPulse.charge for recoPulse in recoPulseList])

            '''
            Removing DOMs with hits less than 200 Hits
            '''
            if sum(recoPulse_chargeList) < self.cuts:
                exit_status = np.array([0])
                exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})
                continue

            '''
            Calculating the mean and removing the tails
            '''

            mean = sum(recoPulse_timeList*recoPulse_chargeList)/sum(recoPulse_chargeList) #mean is weighted
            select_time = recoPulse_timeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
            select_charge = recoPulse_chargeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]

            if len(select_time) < 10 or len(select_charge) < 10:
                exit_status = np.array([1])
                exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})
                continue
            
            exit_status = np.array([3])

            T0 = mean
            D0 = 20.0
            T1 = mean
            D1 = 20.0
            Br = 1.0

            qFunctor_single = LikelihoodFunctor(data,self.pdf,self.time_lim,self.dist_lim)   
          
            minimizer_single = Minuit(qFunctor_single, 
                            t0=T0,
                            error_t0=1.0,
                            limit_t0=(mean-300.,mean+300.),
                            d0=D0,
                            error_d0=1.0,
                            limit_d0=(1.0,120.),
                            t1=T1,
                            error_t1=1.0,
                            limit_t1=(mean-300.,mean+300.),
                            d1=D1,
                            error_d1=1.0,
                            limit_d1=(1.0,120.),
                            br=Br,
                            error_br=1.0,
                            limit_br=(0.0,1.0),
                            errordef=0.5,
                           )

            minimizer_single.fixed["t1"]=True
            minimizer_single.fixed["d1"]=True
            minimizer_single.fixed["br"]=True

            minimizer_single.migrad()

            T0 = mean-20.0
            D0 = 20.0
            T1 = mean+20.0
            D1 = 20.0
            Br = 1.0

            qFunctor_double = LikelihoodFunctor(data,self.pdf,self.time_lim,self.dist_lim)

            minimizer_double = Minuit(qFunctor_double,                                     
                                            t0=T0,
                                            error_t0=1.0,
                                            limit_t0=(mean-300.,mean+300.),
                                            d0=D0,
                                            error_d0=1.0,
                                            limit_d0=(1.0,120.),
                                            t1=T1,
                                            error_t1=1.0,
                                            limit_t1=(mean-300.,mean+300.),
                                            d1=D1,
                                            error_d1=1.0,
                                            limit_d1=(1.0,120.),
                                            br=Br,
                                            error_br=1.0,
                                            limit_br=(0.0,1.0),
                                            errordef=0.5,                                       
                                    )

            minimizer_double.migrad()

            solution_single = minimizer_single.values 
            solution_double = minimizer_double.values

            biGauss_values = np.array([minimizer_single.fval, 
                                       solution_single["t0"],
                                       solution_single["d0"]
                                       ])

            doublePeak_values = np.array([minimizer_double.fval, 
                                          solution_single["t0"],
                                          solution_single["d0"],
                                          solution_single["t1"],
                                          solution_single["d1"],
                                          solution_single["br"]
                                         ])

            biGauss_valuesMap.update({omkey: dataclasses.I3VectorDouble(biGauss_values)})
            doublePeak_valuesMap.update({omkey: dataclasses.I3VectorDouble(doublePeak_values)})
            exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})

        frame[self.output + '_singlePeak'] = biGauss_valuesMap
        frame[self.output + '_doublePeak'] = doublePeak_valuesMap
        frame[self.output + '_exitStatus'] = exitStatusMap
            
        self.PushFrame(frame)    

