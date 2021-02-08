#!/usr/bin/env python                                                                                                
 
# Import some useful ICECUBE modules                                                                                  
from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
import sys
import os
import argparse
import math as m
#import ROOT
from Weights import weight

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")

args = parser.parse_args()

_dir = args.infile 
file_list_aux = os.listdir(_dir)
file_list = [x for x in file_list_aux if '.i3.gz' in x]

nfiles = len(file_list)

#outfile = ROOT.TFile(args.outfile,"RECREATE")
#CorrectTau_DeltaD = ROOT.TH1F("CorrectTau_DeltaD","",500,0.,500.)
#CorrectTau_V1Res = ROOT.TH1F("CorrectTau_V1Res","",500,0.,500.)
#CorrectTau_V2Res = ROOT.TH1F("CorrectTau_V2Res","",500,0.,500.)
#CorrectTau_Energy = ROOT.TH1F("CorrectTau_Energy","",100,0.,10.)
#MissedTau_DeltaD = ROOT.TH1F("MissedTau_DeltaD","",500,0.,500.)
#MissedTau_Energy = ROOT.TH1F("MissedTau_Energy","",100,0.,10.)
#FalseTau_Energy = ROOT.TH1F("FalseTau_Energy","",100.0.,10.)
#CorrectEle_Energy = ROOT.TH1F("CorrectEle_Energy","",100,0,10.)
#CorrectEle_VRes = ROOT.TH1F("CorrectEle_VRes","",500,0.,500.)

#MissIdentified_Tau = ROOT.TH1F("MissIdentified_Tau","",160,-20.,20.)
#CorrIdentified_Tau = ROOT.TH1F("CorrIdentified_Tau","",160,-20.,20.)
#Missed_Tau = ROOT.TH1F("Missed_Tau","",160,-20.,20.)
#CorrectReject = ROOT.TH1F("CorrectReject","",160,-20.,20.)

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

    biGauss_valuesMap = frame['nuTauCurveFit_biGauss']                      
    doublePeak_valuesMap = frame['nuTauCurveFit_doublePeak']

    weightDict = frame["I3MCWeightDict"]
    weight = weightDict['OneWeight']
    #f = astroFlux(NuGPrimary.energy)
    #liveTime = 365*24*60*60
    #event_weight = simpleWeight(weight,f)*liveTime/(2000*weightDict['NEvents'])
    event_weight = weight#(frame, len(file_list))

    mctracktype = MMCTrackList[0].GetI3Particle().type
    #print("Primary Type = "+str(NuGPrimary.type)+" Secondary Type = " 
    #    +str(mctracktype)+" track like = " + str(track_like)+" single like = "
    #    + str(single_like)+" double like = " + str(double_like))

    tau_reco = False
    if double_like < single_like and double_like < track_like :
      tau_reco = True

    likelihood_diff = min(-double_like+single_like,-double_like+track_like)

    maxdif = -100.0
    for key in biGauss_valuesMap.keys() :
      bigauss = biGauss_valuesMap[key]
      doublepeak = doublePeak_valuesMap[key]
      diff = bigauss-doublepeak
      maxdif = max(maxdif,diff)

    true_tau = False
    if NuGPrimary.type == 16 or NuGPrimary.type == -16 :                        
      mctracktype = MMCTrackList[0].GetI3Particle().type      
      if mctracktype == 15 or mctracktype == -15 :
        true_tau = True

    print("Primary Type = "+str(NuGPrimary.type)+" Secondary Type = " +str(mctracktype)+" like dif = "+str(likelihood_diff) + " bigausdiff = "+str(maxdif))  

 #   like = -20.0
 #   while like < likelihood_diff :
 #     if true_tau :
 #       Missed_Tau.Fill(like,event_weight)
 #     else :
 #       CorrectReject.Fill(like,event_weight)
 #     like += 0.25
 #   while like < 20. :
 #     if true_tau :
 #       CorrIdentified_Tau.Fill(like,event_weight)
 #     else :
 #       MissIdentified_Tau.Fill(like,event_weight) 
 #     like += 0.25

#outfile.cd()
#MissIdentified_Tau.Write()           
#CorrIdentified_Tau.Write()               
#Missed_Tau.Write()
#CorrectReject.Write()
#outfile.Close()

    #if NuGPrimary.type == 16 or NuGPrimary.type == -16 :
      #Secondary = MMCTrackList[0].GetI3Particle()
    #  if Secondary.type == 15 or Secondary.type == -15 :
    #    if tau_reco :


      #v1x = MMCTrackList[0].GetXi()
      #v1y = MMCTrackList[0].GetYi()
      #v1z = MMCTrackList[0].GetZi()
      #v2x = MMCTrackList[0].GetXf()
      #v2y = MMCTrackList[0].GetYf()
      #v2z = MMCTrackList[0].GetZf()
      #vertex1 = dataclasses.I3Position(v1x,v1y,v1z)
      #vertex2 = dataclasses.I3Position(v2x,v2y,v2z)

      #energy = MMCTrackList[0].GetEi()
      


