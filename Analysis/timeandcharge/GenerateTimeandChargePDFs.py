from icecube import dataclasses, dataio, simclasses                             
from icecube.icetray import I3Units, OMKey, I3Frame                                    
from icecube.dataclasses import I3Particle
from scipy.stats import norm                                      
import numpy as np                                                              
import sys                                                                      
import os                                                                       
import argparse                                                                 
import math as m 
from ROOT import *

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="",help="")  

args = parser.parse_args()                                                      
                                                                                
_dir = "/data/p-one/tmcelroy/photonbomb/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if '.i3.gz' in x]                         
                                                                                
nfiles = len(file_list)  
print(nfiles)

#file_list = ["/data/p-one/tmcelroy/photonbomb/PhotonBomb_0_prep.i3.gz"]

rootoutfile = TFile(args.outfile,"RECREATE")
TimeChargeDist = TH2F("TimeChargeDist","",400,0.,400.,10010,-10.,10000.)               

PDFBins = []

deltaT = 1.0
deltaD = 1.0

for i in range(400) :
  PDFBins.append([])
  for j in range(10010) :
    PDFBins[-1].append(0.0)

gcd_file = dataio.I3File(args.gcdfile)
geometry = gcd_file.pop_frame()["I3Geometry"]
geoMap = geometry.omgeo  
c = 0.299792458  
n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
c_n = c/n
sigma = 1.0                            
                                                                                
for infile in file_list:                                                        
  infile =  dataio.I3File(os.path.join(_dir,infile))                            
  #infile =  dataio.I3File(infile)
  for frame in infile:                                                          
                                                                                
    if not frame.Has('I3EventHeader') :                                          
      continue

    bomb_pos = frame["PhotonBomb_position"]
    pulse_series = frame["I3Photons"]

    nDOMs = []
    for i in range(400):
      nDOMs.append(0)

    for dom in geoMap.keys() :
      dom_pos = geoMap[dom].position                                      
      distance = np.sqrt((bomb_pos.x-dom_pos.x)**2.0+(bomb_pos.y-dom_pos.y)**2.0+(bomb_pos.z-dom_pos.z)**2.0)
      distbin = int(distance)
      if distbin>-1 and distbin<400 :                                     
        nDOMs[distbin] += 1

    for dom in pulse_series.keys() :
      pulse_list =  pulse_series[dom]
      dom_key = OMKey(dom.string, dom.om, 0)
      dom_pos = geoMap[dom_key].position
      distance = np.sqrt((bomb_pos.x-dom_pos.x)**2.0+(bomb_pos.y-dom_pos.y)**2.0+(bomb_pos.z-dom_pos.z)**2.0)

      for pulse in pulse_list :
        for i in range(-3,4) :
          time = pulse.time-(500.+distance/c_n) + i
          distbin = int(distance) 
          if distbin<0 or distbin>399 :
            continue
          #print("Time = "+str(time) + " distance = "+str(distance)+" norm = "+str(norm(float(i))))
          TimeChargeDist.Fill(float(distance),float(time),TMath.Gaus(float(i),0.0,1.0,True)/float(nDOMs[distbin]))
          timebin = 10+int(time)
          distbin = int(distance)
          if timebin>-1 and timebin<10010 :
            PDFBins[distbin][timebin] += TMath.Gaus(float(i),0.0,1.0,True)/float(nDOMs[distbin])
outfile = file("fittertables.dat","w")
outfile.write("400,10010,0.0,400.0,-10.,10000.\n")

sum_pdf = 0.0
for i in range(400) :                                                           
  for j in range(10010) :
    sum_pdf += PDFBins[i][j]
for i in range(400) :
  for j in range(10010) : 
    outfile.write(str(PDFBins[i][j]/sum_pdf)+"\n")
outfile.close()

rootoutfile.cd()
TimeChargeDist.Write()
rootoutfile.Close()
