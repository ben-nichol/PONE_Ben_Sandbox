from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import I3Particle
from icecube.phys_services import I3Calculator as calc
from icecube.dataclasses import ModuleKey
import numpy as np
import sys, os
import argparse
import math as m
import ROOT


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile",type = str, default="./test_output",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
parser.add_argument("-g", "--gcdfile",type=str,default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")

args = parser.parse_args()

fout = ROOT.TFile(args.outfile+".root","RECREATE")

AttenuationGraphs = []
AttenuationBinCounts = []
for i in range(40) :
    wavelength = 200. + 10.*i
    AttenuationGraphs.append(ROOT.TH1F("AttenuationGraph_"+str(wavelength),"",200,0.0,200))
    AttenuationBinCounts.append(ROOT.TH1F("AttenuationBinCount_"+str(wavelength),"",200,0.0,200))

infilei3 = dataio.I3File(args.infile)

gcdfile = dataio.I3File(args.gcdfile)
for frame in gcdfile :
    if frame.Has('I3Geometry'):
        dom_geo = frame['I3Geometry'].omgeo
        break

for frame in infilei3:

    if not frame.Has('MMCTrackList') :
        continue

    MMCTrackList = frame['MMCTrackList']
    Muon = MMCTrackList[0].GetI3Particle()
    initialpos = [MMCTrackList[0].xi,MMCTrackList[0].yi,MMCTrackList[0].zi]
    finalpos = [MMCTrackList[0].xf,MMCTrackList[0].yf,MMCTrackList[0].zf]
    pulse_series = frame["I3Photons"]
    muon_pos = Muon.pos
    muon_dir = Muon.dir

    for dom in dom_geo.keys() :
        dom_position = dom_geo[dom].position
        reco_dist = calc.cherenkov_distance(Muon, dom_position, 1.35, 1.34)
        clos_app_pos = calc.closest_approach_position(Muon, dom_position)
        l = (clos_app_pos.x-initialpos[0])*muon_dir.x + (clos_app_pos.y-initialpos[1])*muon_dir.y + (clos_app_pos.z-initialpos[2])*muon_dir.z
        if l <= 0 :
            continue
        l = (clos_app_pos.x-finalpos[0])*muon_dir.x + (clos_app_pos.y-finalpos[1])*muon_dir.y + (clos_app_pos.z-finalpos[2])*muon_dir.z
        if l >= 0 :
            continue
        clos_app_dist = np.sqrt((clos_app_pos.x-dom_position.x)**2.0+(clos_app_pos.y-dom_position.y)**2.0+(clos_app_pos.z-dom_position.z)**2.0)
        key = ModuleKey(dom.string,dom.om)
        npulses = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

        if key in pulse_series.keys() :
                if i >= 0 and i < len(npulses) :
                    npulses[i] += max(clos_app_dist,1.0)
                
        for i in range(len(npulses)) :
            AttenuationGraphs[i].Fill(reco_dist,npulses[i])
            AttenuationBinCounts[i].Fill(reco_dist,1.0)

for j in range(len(AttenuationBinCounts)) :
    for i in range(1,201):
        if AttenuationBinCounts[j].GetBinContent(i) > 0.0 :
            AttenuationGraphs[j].SetBinContent(i,AttenuationGraphs[j].GetBinContent(i)/AttenuationBinCounts[j].GetBinContent(i))

fout.cd()
for i in range(len(AttenuationBinCounts)) :
    AttenuationGraphs[i].Write()
fout.Close()
