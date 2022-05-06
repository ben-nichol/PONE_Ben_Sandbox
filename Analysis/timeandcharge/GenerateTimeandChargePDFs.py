from icecube import dataclasses, dataio, simclasses                             
from icecube.icetray import I3Units, OMKey, I3Frame                                    
from icecube.dataclasses import I3Particle
from scipy.stats import norm                                      
import numpy as np                                                              
import sys                                                                      
import os                                                                       
import argparse                                                                 
import math as m 
import pickle

def gaussian(x,x0,sigma):
  return 1./(np.sqrt(2.*np.pi)*sigma)*np.exp(-np.power((x - x0)/sigma, 2.)/2.)

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./output.h5",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="/data/p-one/sim/calibration/sim0012/jobfiles/",help="")  

args = parser.parse_args()                                                      
                                                                                
_dir = args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if '.i3.gz' in x] # and os.stat(os.path.join(_dir,x)).st_size > 800000]                         
                                                                                
nfiles = len(file_list)  
print(nfiles)

dist = list()
char = list()
time = list()
          
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
n = 1.3555714017                                        # 1.33 is the refractive index of water at 20 degrees C
c_n = c/n
sigma = 1.0                           

nDOMs = []
for i in range(400):
  nDOMs.append(0) 
                                                                                
for infile in file_list: 
  print(infile)
  infilei3 =  dataio.I3File(os.path.join(_dir,infile))                            
  for frame in infilei3:                                                          
                                                                                
    if not frame.Has('I3EventHeader') :                                          
      continue

    bomb_pos = frame["PhotonBomb_position"]
    pulse_series = frame["I3Photons"]

    for dom in geoMap.keys() :
      dom_pos = geoMap[dom].position                                      
      distance = np.sqrt((bomb_pos.x-dom_pos.x)**2.0+(bomb_pos.y-dom_pos.y)**2.0+(bomb_pos.z-dom_pos.z)**2.0)
      distbin = int(distance)
      if distbin>-1 and distbin<400 :                                     
        nDOMs[distbin] += 1

    for dom in pulse_series.keys() :
      pulse_list =  pulse_series[dom]
      dom_key = OMKey(dom.string, dom.om, 1)
      dom_pos = geoMap[dom_key].position
      distance = np.sqrt((bomb_pos.x-dom_pos.x)**2.0+(bomb_pos.y-dom_pos.y)**2.0+(bomb_pos.z-dom_pos.z)**2.0)

      for pulse in pulse_list :
        for i in range(-3,4) :
          time = pulse.time-500.+10.0-distance/c_n + i
          distbin = int(distance) 
          if distbin<0 or distbin>399 :
            continue
          timebin = 10+int(time)
          distbin = int(distance)
          if timebin>-1 and timebin<10010 :
            PDFBins[distbin][timebin] += gaussian(float(i),0.0,1.0)

  infilei3.close()

for i in range(400) :
	for j in range(10010) :
		if nDOMs[i] < 1 :
			continue
		PDFBins[i][j] /= nDOMs[i]

for i in range(400) :
        for j in range(10010) :
            dist.append(i)
            char.append(PDFBins[i][j])
            time.append(j)

ouputdict = {"dist":dist,"char":char,"time":time}

pickle.dump(ouputdict,open(args.outfile+".pkl","wb"))

