from icecube import dataclasses, dataio, simclasses                             
from icecube.icetray import I3Units, I3Frame                                    
from icecube.dataclasses import I3Particle
from scipy.stats import norm                                      
import numpy as np                                                              
import sys                                                                      
import os                                                                       
import argparse                                                                 
import math as m 

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
                                                       
args = parser.parse_args()                                                      
                                                                                
_dir = args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if '.i3.gz' in x]                         
                                                                                
nfiles = len(file_list)  

outfile = ROOT.TFile(args.outfile,"RECREATE")
TimeChargeDist = ROOT.TH1F("TimeChargeDist","",200,0.,200.,2000,-500.,1500.)               

gcd_file = dataio.I3File(args.gcdfile)
geometry = gcd_file.pop_frame()["I3Geometry"]
geoMap = geometry.omgeo  
c = 0.299792458  
n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
c_n = c/n
sigma = 1.0                            
                                                                                
for infile in file_list:                                                        
  infile =  dataio.I3File(os.path.join(_dir,infile))                            
  for frame in infile:                                                          
                                                                                
    if not frame.Has('MCEventHeader') :                                          
      continue

    bomb_pos = frame["PhotonBomb_position"]
    pulse_series = frame["I3Photons"]

    for dom in pulse_series :
    	pulse_list =  pulse_series[dom]
    	dom_pos = geoMap[dom].position
    	distance = np.sqrt((bomb_pos.x-dom_pos.x)**2.0+(bomb_pos.y-dom_pos.y)**2.0+(bomb_pos.z-dom_pos.z)**2.0)
    	for pulse in pulse_list :
        for i in range(-3,4) :
    		  TimeChargeDist.Fill(distance,pulse.time-(500.+distance/c_n)*I3Units.ns,norm(i))

outfile = file("fittertables.dat","w")
outfile.write("200,2000,0.0,200.0,-500.,1500.\n")
for i in range(TimeChargeDist.GetN()) :
	outfile.Write(str(TimeChargeDist.GetBinContent(1+i))+"\n")
outfile.close()

outfile.cd()
TimeChargeDist.write()
outfile.Close()