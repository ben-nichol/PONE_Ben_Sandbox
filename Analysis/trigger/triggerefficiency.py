from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame 
from icecube.dataclasses import I3Particle 
import numpy as np                 
import sys, os
import argparse
import math as m
from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="",help="")  

args = parser.parse_args()    

_dir = "/data/p-one/tmcelroy/muons/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if ('.i3.gz' in x and 'TrigReco_Corsika__PhotonProp' in x)]      
gcd_file = dataio.I3File("/home/users/tmcelroy/pone_offline/GCD/PONE_Phase1.i3.gz")
for frame in gcd_file:
    domsUsed = frame['I3Geometry'].omgeo

DetectorTrigger_Energy = list()
DetectorTrigger_Zenith = list()
StringTrigger_Energy = list()
StringTrigger_Zenith = list()
InterStringTrigger_Energy = list()
InterStringTrigger_Zenith = list()
AllEvents_Energy = list()
AllEvents_Zenith = list()

filecount = 0

for infile in file_list :
    print(infile)
    infilei3 = dataio.I3File(os.path.join(_dir,infile))  
    filecount += 1
    if filecount> 20 :
        break
    for frame in infilei3:                                                                                                                                      
        if not frame.Has('I3EventHeader') :                                          
            continue

        detectorTriggerTime = frame["DetectorTriggers"]
        stringTriggerTime = frame["SingleStringTriggers"]
        interStringTriggerTime = frame["InterStringTriggers"]

        MMCTrackList = frame['MMCTrackList']
        Muon = MMCTrackList[0].GetI3Particle()
        #pulse_series = frame["SignificanHits"]

        AllEvents_Energy.append(Muon.energy)
        AllEvents_Zenith.append(Muon.dir.zenith)

        if len(detectorTriggerTime) > 0 :
            DetectorTrigger_Energy.append(Muon.energy)
            DetectorTrigger_Zenith.append(Muon.dir.zenith)
        if len(stringTriggerTime) > 0 :
            StringTrigger_Energy.append(Muon.energy)
            StringTrigger_Zenith.append(Muon.dir.zenith)
        if len(interStringTriggerTime) > 0 :
            InterStringTrigger_Energy.append(Muon.energy)
            InterStringTrigger_Zenith.append(Muon.dir.zenith)

np_DetectorTrigger_Energy = np.array(DetectorTrigger_Energy)
np_DetectorTrigger_Zenith = np.array(DetectorTrigger_Zenith)
np_StringTrigger_Energy = np.array(StringTrigger_Energy)
np_StringTrigger_Zenith = np.array(StringTrigger_Zenith)
np_InterStringTrigger_Energy = np.array(InterStringTrigger_Energy)
np_InterStringTrigger_Zenith = np.array(InterStringTrigger_Zenith)
np_AllEvents_Energy = np.array(AllEvents_Energy)
np_AllEvents_Zenith = np.array(AllEvents_Zenith)

n, bins, patches = plt.hist([np_AllEvents_Energy,np_DetectorTrigger_Energy], 
                            50 ,
                            range = (min(AllEvents_Energy),max(AllEvents_Energy)),
                            label=["All Events","3 DOM Coincidence"])

plt.xlabel('Energy (GeV)')
plt.ylabel('N Triggers')
plt.title('3 DOM Coincidence')
plt.grid(True)
plt.savefig("Detector_Energy.png", format="png")
plt.clf() 

eff = list()
E = list()
for i in range(len(n[0])) :
    if n[0][i] > 0 :
        eff.append(n[1][i]/n[0][i])
    E.append(0.5*(bins[0]+bins[1]))
print(eff)
plt.plot(E,eff)
plt.xlabel('Energy (GeV)')
plt.ylabel('Trigger Eff')
plt.savefig("Detector_Energy_Eff.png", format="png")
plt.clf()

n, bins, patches = plt.hist([np_AllEvents_Zenith,np_DetectorTrigger_Zenith], 
                            50 ,
                            range = (min(AllEvents_Zenith),max(AllEvents_Zenith)),
                            label=["All Events","3 DOM Coincidence"])

plt.xlabel('Zenith')
plt.ylabel('N Triggers')
plt.title('3 DOM Coincidence')
plt.grid(True)
plt.savefig("Detector_Zenith.png", format="png")
plt.clf()

n, bins, patches = plt.hist([np_AllEvents_Energy,np_StringTrigger_Energy], 
                            50 ,
                            range = (min(AllEvents_Energy),max(AllEvents_Energy)),
                            label=["All Events","2 Adjacent String DOMs"])

plt.xlabel('Energy (GeV)')
plt.ylabel('N Triggers')
plt.title('2 Adjacent String DOMs')
plt.grid(True)
plt.savefig("SingleString_Energy.png", format="png")
plt.clf()

n, bins, patches = plt.hist([np_AllEvents_Zenith,np_StringTrigger_Zenith],
                            50 ,
                            range = (min(AllEvents_Zenith),max(AllEvents_Zenith)),
                            label=["All Events","2 Adjacent String DOMs"])

plt.xlabel('Zenith')
plt.ylabel('N Triggers')
plt.title('2 Adjacent String DOMs')
plt.grid(True)
plt.savefig("SingleString_Zenith.png", format="png")
plt.clf()

n, bins, patches = plt.hist([np_AllEvents_Energy,np_InterStringTrigger_Energy],
                            50 ,
                            range = (min(AllEvents_Energy),max(AllEvents_Energy)),
                            label=["All Events","2 Adjacent Same Z DOMs"])

plt.xlabel('Energy (GeV)')
plt.ylabel('N Triggers')
plt.title('2 Adjacent Same Z DOMs')
plt.grid(True)
plt.savefig("InterString_Energy.png", format="png")
plt.clf()

n, bins, patches = plt.hist([np_AllEvents_Zenith,np_InterStringTrigger_Zenith],
                            50 ,
                            range = (min(AllEvents_Zenith),max(AllEvents_Zenith)),
                            label=["All Events","2 Adjacent Same Z DOMs"])

plt.xlabel('Zenith')
plt.ylabel('N Triggers')
plt.title('2 Adjacent Same Z DOMs')
plt.grid(True)
plt.savefig("InterString_Zenith.png", format="png")
plt.clf()

#fout.cd()
#fout.Close()










