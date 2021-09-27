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
import ROOT 

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="",help="")  

args = parser.parse_args()    

_dir = "/data/p-one/tmcelroy/muons/TriggerTest/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if ('.i3.gz' in x and 'TrigReco_Corsika_lowE_PhotonProp_' in x)]      
gcd_file = dataio.I3File("/home/users/tmcelroy/pone_offline/GCD/PONE_Phase1.i3.gz")
for frame in gcd_file:
    domsUsed = frame['I3Geometry'].omgeo

filecount = 0

fout = ROOT.TFile.Open("TriggerEfficiency.root","RECREATE")

nbinenergy = 240
nbinzenith = 200
nbinsazimuth = 200
nbinscharge = 100
nbinsdoms = 100
nbinscoinc = 20

minenergy = 0.0
maxenergy = 24.0

minzenith = -1.0
maxzenith = 1.0

minazimuth = 0.0
maxazimuth = 2.0*np.pi

mincharge = 0.0
maxcharge = 1000

minndoms = 0
maxndoms = 100

mincoinc = 0
maxcoinc = 20

AllEvents_Energy_Hist = ROOT.TH1F("AllEvents_Energy_Hist","",nbinenergy,minenergy,maxenergy)
AllEvents_Zenith_Hist = ROOT.TH1F("AllEvents_Zenith_Hist","",nbinzenith,minzenith,maxzenith)
AllEvents_Azimuth_Hist = ROOT.TH1F("AllEvents_Azimuth_Hist","",nbinsazimuth,minazimuth,maxazimuth)

DetectorTrigger_Energy_3PMT_2DOM_Hist = ROOT.TH1F("DetectorTrigger_Energy_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy)
DetectorTrigger_Zenith_3PMT_2DOM_Hist = ROOT.TH1F("DetectorTrigger_Zenith_3PMT_2DOM_Hist","",nbinzenith,minzenith,maxzenith)
DetectorTrigger_Azimuth_3PMT_2DOM_Hist = ROOT.TH1F("DetectorTrigger_Azimuth_3PMT_2DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
StringTrigger_Energy_3PMT_2DOM_Hist = ROOT.TH1F("StringTrigger_Energy_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy)
StringTrigger_Zenith_3PMT_2DOM_Hist = ROOT.TH1F("StringTrigger_Zenith_3PMT_2DOM_Hist","",nbinzenith,minzenith,maxzenith)
StringTrigger_Azimuth_3PMT_2DOM_Hist = ROOT.TH1F("StringTrigger_Azimuth_3PMT_2DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
AllTriggers_Energy_3PMT_2DOM_Hist = ROOT.TH1F("AllTriggers_Energy_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy)
AllTriggers_Zenith_3PMT_2DOM_Hist = ROOT.TH1F("AllTriggers_Zenith_3PMT_2DOM_Hist","",nbinzenith,minzenith,maxzenith)
AllTriggers_Azimuth_3PMT_2DOM_Hist = ROOT.TH1F("AllTriggers_Azimuth_3PMT_2DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)

MissingTriggerCharge_3PMT_2DOM_Hist = ROOT.TH2F("MissingTriggerCharge_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscharge,mincharge,maxcharge)
MissingTriggerCharge_2PMT_4DOM_Hist = ROOT.TH2F("MissingTriggerCharge_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscharge,mincharge,maxcharge)
MissingTriggerCharge_2PMT_3DOM_Hist = ROOT.TH2F("MissingTriggerCharge_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscharge,mincharge,maxcharge)

MissingTriggerNDOMs_3PMT_2DOM_Hist = ROOT.TH2F("MissingTriggerNDOMs_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinsdoms,minndoms,maxndoms)
MissingTriggerNDOMs_2PMT_4DOM_Hist = ROOT.TH2F("MissingTriggerNDOMs_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinsdoms,minndoms,maxndoms)
MissingTriggerNDOMs_2PMT_3DOM_Hist = ROOT.TH2F("MissingTriggerNDOMs_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinsdoms,minndoms,maxndoms)

MissingTriggerNCoinc_3PMT_2DOM_Hist = ROOT.TH2F("MissingTriggerNCoinc_3PMT_2DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscoinc,mincoinc,maxcoinc)
MissingTriggerNCoinc_2PMT_4DOM_Hist = ROOT.TH2F("MissingTriggerNCoinc_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscoinc,mincoinc,maxcoinc)
MissingTriggerNCoinc_2PMT_3DOM_Hist = ROOT.TH2F("MissingTriggerNCoinc_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy,nbinscoinc,mincoinc,maxcoinc)

DetectorTrigger_Energy_2PMT_4DOM_Hist = ROOT.TH1F("DetectorTrigger_Energy_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy)
DetectorTrigger_Zenith_2PMT_4DOM_Hist = ROOT.TH1F("DetectorTrigger_Zenith_2PMT_4DOM_Hist","",nbinzenith,minzenith,maxzenith)
DetectorTrigger_Azimuth_2PMT_4DOM_Hist = ROOT.TH1F("DetectorTrigger_Azimuth_2PMT_4DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
StringTrigger_Energy_2PMT_4DOM_Hist = ROOT.TH1F("StringTrigger_Energy_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy)
StringTrigger_Zenith_2PMT_4DOM_Hist = ROOT.TH1F("StringTrigger_Zenith_2PMT_4DOM_Hist","",nbinzenith,minzenith,maxzenith)
StringTrigger_Azimuth_2PMT_4DOM_Hist = ROOT.TH1F("StringTrigger_Azimuth_2PMT_4DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
AllTriggers_Energy_2PMT_4DOM_Hist = ROOT.TH1F("AllTriggers_Energy_2PMT_4DOM_Hist","",nbinenergy,minenergy,maxenergy)
AllTriggers_Zenith_2PMT_4DOM_Hist = ROOT.TH1F("AllTriggers_Zenith_2PMT_4DOM_Hist","",nbinzenith,minzenith,maxzenith)
AllTriggers_Azimuth_2PMT_4DOM_Hist = ROOT.TH1F("AllTriggers_Azimuth_2PMT_4DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)

DetectorTrigger_Energy_2PMT_3DOM_Hist = ROOT.TH1F("DetectorTrigger_Energy_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy)
DetectorTrigger_Zenith_2PMT_3DOM_Hist = ROOT.TH1F("DetectorTrigger_Zenith_2PMT_3DOM_Hist","",nbinzenith,minzenith,maxzenith)
DetectorTrigger_Azimuth_2PMT_3DOM_Hist = ROOT.TH1F("DetectorTrigger_Azimuth_2PMT_3DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
StringTrigger_Energy_2PMT_3DOM_Hist = ROOT.TH1F("StringTrigger_Energy_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy)
StringTrigger_Zenith_2PMT_3DOM_Hist = ROOT.TH1F("StringTrigger_Zenith_2PMT_3DOM_Hist","",nbinzenith,minzenith,maxzenith)
StringTrigger_Azimuth_2PMT_3DOM_Hist = ROOT.TH1F("StringTrigger_Azimuth_2PMT_3DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)
AllTriggers_Energy_2PMT_3DOM_Hist = ROOT.TH1F("AllTriggers_Energy_2PMT_3DOM_Hist","",nbinenergy,minenergy,maxenergy)
AllTriggers_Zenith_2PMT_3DOM_Hist = ROOT.TH1F("AllTriggers_Zenith_2PMT_3DOM_Hist","",nbinzenith,minzenith,maxzenith)
AllTriggers_Azimuth_2PMT_3DOM_Hist = ROOT.TH1F("AllTriggers_Azimuth_2PMT_3DOM_Hist","",nbinsazimuth,minazimuth,maxazimuth)

def indetector(d,p):
    R_cyl = 150.
    Z_cyl = 400.

    #Never hit within cylinder
    l = -(p.x*d.x+p.y*d.y)/(d.x**2.0+d.y**2.0)
    r2 = (p.x+l*d.x)**2.0+(p.y+l*d.y)**2.0
    if r2 > R_cyl**2.0 :
        return False

    #The closest point is within z limits.
    z_closest = p.z+l*d.z
    if z_closest < Z_cyl and z_closest > -Z_cyl :
        return True
    else: False

    #When it hist the top and bottom, is it within radius?
    l_upper = (Z_cyl-p.z)/d.z
    r2_upper = (p.x+l_upper*d.x)**2.0+(p.y+l_upper*d.y)**2.0
    l_lower = (-Z_cyl-p.z)/d.z
    r2_lower = (p.x+l_lower*d.x)**2.0+(p.y+l_lower*d.y)**2.0

    if r2_upper > R_cyl**2.0 or r2_lower > R_cyl**2.0 :
        return False

    return True


for infile in file_list :
    print(infile)
    infilei3 = dataio.I3File(os.path.join(_dir,infile))  
    for frame in infilei3:                                                                                                                                      
        if not frame.Has('I3EventHeader') :                                          
            continue

        detectorTriggerTime_3PMT_2DOM = frame["DetectorTriggers_3PMT_2DOM"]
        stringTriggerTime_3PMT_2DOM = frame["StringTriggers_3PMT_2DOM"]

        detectorTriggerTime_2PMT_4DOM = frame["DetectorTriggers_2PMT_4DOM"]
        stringTriggerTime_2PMT_4DOM = frame["StringTriggers_2PMT_4DOM"]

        detectorTriggerTime_2PMT_3DOM = frame["DetectorTriggers_2PMT_3DOM"]
        stringTriggerTime_2PMT_3DOM = frame["StringTriggers_2PMT_3DOM"]


        MMCTrackList = frame['MMCTrackList']
        Muon = MMCTrackList[0].GetI3Particle()
        
        if not indetector(Muon.dir, Muon.pos) :
            continue

        totalcharge = 0.0
        DOMList = {}
        pulseseries = frame["I3Photons_PMTResponse"]

        for key in pulseseries.keys() :
            omkey = OMKey(key.string,key.om,0)
            if omkey not in DOMList.keys() :
                DOMList[omkey] = 1
            else :
                DOMList[omkey] += 1
            for pulse in pulseseries[key] :
                totalcharge += pulse.charge

        if totalcharge < 20.0 :
            continue
        if len(DOMList.keys()) < 3.0 :
            continue

        maxcoinc = 0

        for dom in DOMList :
            if DOMList[dom] > maxcoinc :
                maxcoinc = DOMList[dom]

        #pulse_series = frame["SignificanHits"]

        logE = ROOT.TMath.Log10(Muon.energy)
        azimuth = Muon.dir.azimuth
        cos_zenith = ROOT.TMath.Cos(Muon.dir.zenith)

        AllEvents_Energy_Hist.Fill(logE)
        AllEvents_Zenith_Hist.Fill(cos_zenith)
        AllEvents_Azimuth_Hist.Fill(azimuth) 

        #print(ROOT.TMath.Log(Muon.energy))

        ntrigger_3PMT_2DOM = 0
        if len(detectorTriggerTime_3PMT_2DOM) > 0 :
            DetectorTrigger_Energy_3PMT_2DOM_Hist.Fill(logE)
            DetectorTrigger_Zenith_3PMT_2DOM_Hist.Fill(cos_zenith)
            DetectorTrigger_Azimuth_3PMT_2DOM_Hist.Fill(azimuth)
            ntrigger_3PMT_2DOM += 1
        if len(stringTriggerTime_3PMT_2DOM) > 0 :
            StringTrigger_Energy_3PMT_2DOM_Hist.Fill(logE)
            StringTrigger_Zenith_3PMT_2DOM_Hist.Fill(cos_zenith)
            StringTrigger_Azimuth_3PMT_2DOM_Hist.Fill(azimuth)
            ntrigger_3PMT_2DOM += 1
        if ntrigger_3PMT_2DOM > 0 :
            AllTriggers_Energy_3PMT_2DOM_Hist.Fill(logE)
            AllTriggers_Zenith_3PMT_2DOM_Hist.Fill(cos_zenith)
            AllTriggers_Azimuth_3PMT_2DOM_Hist.Fill(azimuth)

        ntrigger_2PMT_4DOM = 0
        if len(detectorTriggerTime_2PMT_4DOM) > 0 :
            DetectorTrigger_Energy_2PMT_4DOM_Hist.Fill(logE)
            DetectorTrigger_Zenith_2PMT_4DOM_Hist.Fill(cos_zenith)
            DetectorTrigger_Azimuth_2PMT_4DOM_Hist.Fill(azimuth)
            ntrigger_2PMT_4DOM += 1
        if len(stringTriggerTime_2PMT_4DOM) > 0 :
            StringTrigger_Energy_2PMT_4DOM_Hist.Fill(logE)
            StringTrigger_Zenith_2PMT_4DOM_Hist.Fill(cos_zenith)
            StringTrigger_Azimuth_2PMT_4DOM_Hist.Fill(azimuth)
            ntrigger_2PMT_4DOM += 1
        if ntrigger_2PMT_4DOM > 0 :
            AllTriggers_Energy_2PMT_4DOM_Hist.Fill(logE)
            AllTriggers_Zenith_2PMT_4DOM_Hist.Fill(cos_zenith)
            AllTriggers_Azimuth_2PMT_4DOM_Hist.Fill(azimuth)

        ntrigger_2PMT_3DOM = 0
        if len(detectorTriggerTime_2PMT_3DOM) > 0 :
            DetectorTrigger_Energy_2PMT_3DOM_Hist.Fill(logE)
            DetectorTrigger_Zenith_2PMT_3DOM_Hist.Fill(cos_zenith)
            DetectorTrigger_Azimuth_2PMT_3DOM_Hist.Fill(azimuth)
            ntrigger_2PMT_3DOM += 1
        if len(stringTriggerTime_2PMT_3DOM) > 0 :
            StringTrigger_Energy_2PMT_3DOM_Hist.Fill(logE)
            StringTrigger_Zenith_2PMT_3DOM_Hist.Fill(cos_zenith)
            StringTrigger_Azimuth_2PMT_3DOM_Hist.Fill(azimuth)
            ntrigger_2PMT_3DOM += 1
        if ntrigger_2PMT_3DOM > 0 :
            AllTriggers_Energy_2PMT_3DOM_Hist.Fill(logE)
            AllTriggers_Zenith_2PMT_3DOM_Hist.Fill(cos_zenith)
            AllTriggers_Azimuth_2PMT_3DOM_Hist.Fill(azimuth)

        if (ntrigger_3PMT_2DOM == 0) or (ntrigger_2PMT_4DOM == 0) or (ntrigger_2PMT_3DOM == 0) :
            
            if ntrigger_3PMT_2DOM == 0 :
                MissingTriggerCharge_3PMT_2DOM_Hist.Fill(logE,totalcharge)
                MissingTriggerNDOMs_3PMT_2DOM_Hist.Fill(logE,len(DOMList))
                MissingTriggerNCoinc_3PMT_2DOM_Hist.Fill(logE,maxcoinc)
            if ntrigger_2PMT_4DOM == 0 :
                MissingTriggerCharge_2PMT_4DOM_Hist.Fill(logE,totalcharge)
                MissingTriggerNDOMs_2PMT_4DOM_Hist.Fill(logE,len(DOMList))
                MissingTriggerNCoinc_2PMT_4DOM_Hist.Fill(logE,maxcoinc)
            if ntrigger_2PMT_3DOM == 0 :
                MissingTriggerCharge_2PMT_3DOM_Hist.Fill(logE,totalcharge)
                MissingTriggerNDOMs_2PMT_3DOM_Hist.Fill(logE,len(DOMList))
                MissingTriggerNCoinc_2PMT_3DOM_Hist.Fill(logE,maxcoinc)


for i in range(1,nbinenergy+1):
    if AllEvents_Energy_Hist.GetBinContent(i) > 0.0 :
        DetectorTrigger_Energy_3PMT_2DOM_Hist.SetBinContent(i,DetectorTrigger_Energy_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        StringTrigger_Energy_3PMT_2DOM_Hist.SetBinContent(i,StringTrigger_Energy_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        AllTriggers_Energy_3PMT_2DOM_Hist.SetBinContent(i,AllTriggers_Energy_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))

        DetectorTrigger_Energy_2PMT_4DOM_Hist.SetBinContent(i,DetectorTrigger_Energy_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        StringTrigger_Energy_2PMT_4DOM_Hist.SetBinContent(i,StringTrigger_Energy_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        AllTriggers_Energy_2PMT_4DOM_Hist.SetBinContent(i,AllTriggers_Energy_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))

        DetectorTrigger_Energy_2PMT_3DOM_Hist.SetBinContent(i,DetectorTrigger_Energy_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        StringTrigger_Energy_2PMT_3DOM_Hist.SetBinContent(i,StringTrigger_Energy_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))
        AllTriggers_Energy_2PMT_3DOM_Hist.SetBinContent(i,AllTriggers_Energy_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Energy_Hist.GetBinContent(i))

for i in range(1,nbinzenith+1):
    if AllEvents_Zenith_Hist.GetBinContent(i) > 0.0 :
        DetectorTrigger_Zenith_3PMT_2DOM_Hist.SetBinContent(i,DetectorTrigger_Zenith_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        StringTrigger_Zenith_3PMT_2DOM_Hist.SetBinContent(i,StringTrigger_Zenith_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        AllTriggers_Zenith_3PMT_2DOM_Hist.SetBinContent(i,AllTriggers_Zenith_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))

        DetectorTrigger_Zenith_2PMT_4DOM_Hist.SetBinContent(i,DetectorTrigger_Zenith_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        StringTrigger_Zenith_2PMT_4DOM_Hist.SetBinContent(i,StringTrigger_Zenith_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        AllTriggers_Zenith_2PMT_4DOM_Hist.SetBinContent(i,AllTriggers_Zenith_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))

        DetectorTrigger_Zenith_2PMT_3DOM_Hist.SetBinContent(i,DetectorTrigger_Zenith_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        StringTrigger_Zenith_2PMT_3DOM_Hist.SetBinContent(i,StringTrigger_Zenith_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))
        AllTriggers_Zenith_2PMT_3DOM_Hist.SetBinContent(i,AllTriggers_Zenith_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Zenith_Hist.GetBinContent(i))

for i in range(1,nbinsazimuth+1):                                                                                                                                         
    if AllEvents_Azimuth_Hist.GetBinContent(i) > 0.0 : 
        DetectorTrigger_Azimuth_3PMT_2DOM_Hist.SetBinContent(i,DetectorTrigger_Azimuth_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        StringTrigger_Azimuth_3PMT_2DOM_Hist.SetBinContent(i,StringTrigger_Azimuth_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        AllTriggers_Azimuth_3PMT_2DOM_Hist.SetBinContent(i,AllTriggers_Azimuth_3PMT_2DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))                         

        DetectorTrigger_Azimuth_2PMT_4DOM_Hist.SetBinContent(i,DetectorTrigger_Azimuth_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        StringTrigger_Azimuth_2PMT_4DOM_Hist.SetBinContent(i,StringTrigger_Azimuth_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        AllTriggers_Azimuth_2PMT_4DOM_Hist.SetBinContent(i,AllTriggers_Azimuth_2PMT_4DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))

        DetectorTrigger_Azimuth_2PMT_3DOM_Hist.SetBinContent(i,DetectorTrigger_Azimuth_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        StringTrigger_Azimuth_2PMT_3DOM_Hist.SetBinContent(i,StringTrigger_Azimuth_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))
        AllTriggers_Azimuth_2PMT_3DOM_Hist.SetBinContent(i,AllTriggers_Azimuth_2PMT_3DOM_Hist.GetBinContent(i)/AllEvents_Azimuth_Hist.GetBinContent(i))

fout.cd()

DetectorTrigger_Energy_3PMT_2DOM_Hist.Write()
DetectorTrigger_Zenith_3PMT_2DOM_Hist.Write()
DetectorTrigger_Azimuth_3PMT_2DOM_Hist.Write()
StringTrigger_Energy_3PMT_2DOM_Hist.Write()
StringTrigger_Zenith_3PMT_2DOM_Hist.Write()
StringTrigger_Azimuth_3PMT_2DOM_Hist.Write()
AllTriggers_Energy_3PMT_2DOM_Hist.Write()
AllTriggers_Zenith_3PMT_2DOM_Hist.Write()
AllTriggers_Azimuth_3PMT_2DOM_Hist.Write()

DetectorTrigger_Energy_2PMT_4DOM_Hist.Write()
DetectorTrigger_Zenith_2PMT_4DOM_Hist.Write()
DetectorTrigger_Azimuth_2PMT_4DOM_Hist.Write()
StringTrigger_Energy_2PMT_4DOM_Hist.Write()
StringTrigger_Zenith_2PMT_4DOM_Hist.Write()
StringTrigger_Azimuth_2PMT_4DOM_Hist.Write()
AllTriggers_Energy_2PMT_4DOM_Hist.Write()
AllTriggers_Zenith_2PMT_4DOM_Hist.Write()
AllTriggers_Azimuth_2PMT_4DOM_Hist.Write()

DetectorTrigger_Energy_2PMT_3DOM_Hist.Write()
DetectorTrigger_Zenith_2PMT_3DOM_Hist.Write()
DetectorTrigger_Azimuth_2PMT_3DOM_Hist.Write()
StringTrigger_Energy_2PMT_3DOM_Hist.Write()
StringTrigger_Zenith_2PMT_3DOM_Hist.Write()
StringTrigger_Azimuth_2PMT_3DOM_Hist.Write()
AllTriggers_Energy_2PMT_3DOM_Hist.Write()
AllTriggers_Zenith_2PMT_3DOM_Hist.Write()
AllTriggers_Azimuth_2PMT_3DOM_Hist.Write()

MissingTriggerCharge_3PMT_2DOM_Hist.Write()
MissingTriggerNDOMs_3PMT_2DOM_Hist.Write()
MissingTriggerCharge_2PMT_4DOM_Hist.Write()
MissingTriggerNDOMs_2PMT_4DOM_Hist.Write()
MissingTriggerCharge_2PMT_3DOM_Hist.Write()
MissingTriggerNDOMs_2PMT_3DOM_Hist.Write()
MissingTriggerNCoinc_3PMT_2DOM_Hist.Write()
MissingTriggerNCoinc_2PMT_4DOM_Hist.Write()
MissingTriggerNCoinc_2PMT_3DOM_Hist.Write()

AllEvents_Energy_Hist.Write()
AllEvents_Zenith_Hist.Write()
AllEvents_Azimuth_Hist.Write()

fout.Close()


