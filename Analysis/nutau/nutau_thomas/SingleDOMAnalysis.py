from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
import sys
import os
import argparse
import math as m
import ROOT
from Weights import weight

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")

args = parser.parse_args()

_dir = args.infile 
file_list_aux = os.listdir(_dir)
file_list = [x for x in file_list_aux if '.i3.gz' in x]

nfiles = len(file_list)

outfile = ROOT.TFile(args.outfile,"RECREATE")

MissIdentified_Tau = ROOT.TH1F("MissIdentified_Tau","",40,0.,20.)
CorrIdentified_Tau = ROOT.TH1F("CorrIdentified_Tau","",40,0.,20.)
Missed_Tau = ROOT.TH1F("Missed_Tau","",40,0.,20.)

for infile in file_list:
  infile =  dataio.I3File(os.path.join(_dir,infile))  
  for frame in infile:

    if not frame.Has('MMCTrackList') :
      continue

    MMCTrackList = frame['MMCTrackList']
    if len(MMCTrackList)<1 :
      continue
    NuGPrimary = frame['NuGPrimary']

    track_like = frame["llhfit_nloglike"].value
    single_like = frame["NuTau_single_nlogl"].value
    double_like = frame["NuTau_double_nlogl"].value

    biGauss_valuesMap = frame[self.output+'_biGauss']
    doublePeak_valuesMap = frame[self.output+ '_doublePeak']

    for dom in biGauss_valuesMap :
    	diff = biGauss_valuesMap[dom] - doublePeak_valuesMap[dom]

    event_weight = weight(frame,nfiles)

    mctracktype = MMCTrackList[0].GetI3Particle().type
    print("Primary Type = "+str(NuGPrimary.type)+" Secondary Type = " 
        +str(mctracktype)+" track like = " + str(track_like)+" single like = "
        + str(single_like)+" double like = " + str(double_like))

    tau_reco = False
    if double_like < single_like and double_like < track_like :
      tau_reco = True

    likelihood_diff = min(double_like-single_like,double_like-track_like)

    true_tau = False
    if NuGPrimary.type == 16 or NuGPrimary.type == -16 :                        
      mctracktype = MMCTrackList[0].GetI3Particle().type      
      if mctracktype == 15 or mctracktype == -15 :
        true_tau = True

    like = 0.0
    while like < likelihood_diff :
      if true_tau :
        Missed_Tau.Fill(like,event_weight)
      else :
        MissIdentified_Tau.Fill(like,event_weight)
      like += 0.5
    while like < 20. :
      if true_tau :
        CorrIdentified_Tau.Fill(like,event_weight)
      like += 0.5

outfile.cd()
MissIdentified_Tau.Write()           
CorrIdentified_Tau.Write()               
Missed_Tau.Write()
outfile.Close()
