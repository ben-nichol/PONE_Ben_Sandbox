from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np                 
import sys
import argparse
import math as m
import ROOT

parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="",help="")  

args = parser.parse_args()    

fout = ROOT.TFile("MuonRecoRes.root","RECREATE")  
DirectionRes_linefit = TH1F("DirectionRes_linefit","",100,-1.0,1.0)
DirectionRes_llhfit = TH1F("DirectionRes_llhfit","",100,-1.0,1.0)
VertexRes_llhfit = TH1F("VertexRes_llhfit","",100,0.0,50.0)
VertexRes_linefit = TH1F("VertexRes_linefit","",100,0.0,50.0)
DirectionRes_linefit_vs_minr = TH2F("DirectionRes_linefit_vs_minr","",110,0.,550.,100,-1.0,1.0)
DirectionRes_llhfit_vs_minr = TH2F("DirectionRes_llhfit_vs_minr","",110,0.,550.,100,-1.0,1.0)

_dir = "/data/p-one/tmcelroy/muons/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if ('.i3.gz' in x and 'NuTauFit_Corsika__PhotonProp' in x)]      

for infile in file_list :
	infilei3 = dataio.I3File(os.path.join(_dir,infile))  

	for frame in infile:                                                          
                                                                                
    	if not frame.Has('I3EventHeader') :                                          
      		continue

      	MMCTrackList = frame['MMCTrackList']
      	Muon = MMCTrackList[0].GetI3Particle()

      	p_2 = Muon.pos.x**2+Muon.pos.y**2+Muon.pos.z**2
      	pd = -(Muon.pos.x*Muon.dir.x+Muon.pos.y*Muon.dir.y+Muon.pos.z*Muon.dir.z)
      	r_2 = 550.0**2

      	L = -2.0*pd - np.sqrt(pd**2.0-p_2+r_2)

      	muon_vertex = dataclasses.I3Position(Muon.pos.x-L*Muon.dir.x,Muon.pos.y-L*Muon.dir.y,Muon.pos.z-L*Muon.dir.z)
      	muon_time = Muon.time -L/0.299792458
      	muon_direction = dataclasses.I3Direction(-Muon.dir.x,-Muon.dir.y,-Muon.dir.z)

      	#Closest Approach
      	minl = -pd
      	minpos = dataclasses.I3Position(Muon.pos.x-minl*Muon.dir.x,Muon.pos.y-minl*Muon.dir.y,Muon.pos.z-minl*Muon.dir.z)
      	minr = np.sqrt(minpos.x**2+minpos.y**2+minpos.z**2)

      	#within detector
      	min_cyl_r = np.sqrt(minpos.x**2+minpos.y**2)
      	min_z = minpos.z

      	if np.abs(min_z) > 500. or min_cyl_r > 200. :
      		continue

      	linefit = frame['linefit']
      	llhfit = frame['llhfit']

      	solidangle_linefit = muon_direction.x*linefit.dir.x+muon_direction.y*linefit.dir.y+muon_direction.z*linefit.dir.z
      	solidangle_llhfit = muon_direction.x*llhfit.dir.x+muon_direction.y*llhfit.dir.y+muon_direction.z*llhfit.dir.z
      	vertexdiff_linefit = np.sqrt((muon_vertex.x-linefit.pos.x)**2.0+(muon_vertex.y-linefit.pos.y)**2.0+(muon_vertex.z-linefit.pos.z)**2.0)
      	vertexdiff_llhfit = np.sqrt((muon_vertex.x-llhfit.pos.x)**2.0+(muon_vertex.y-llhfit.pos.y)**2.0+(muon_vertex.z-llhfit.pos.z)**2.0)

      	DirectionRes_linefit.Fill(solidangle_linefit)
      	DirectionRes_llhfit.Fill(solidangle_llhfit)
		VertexRes_llhfit.Fill(vertexdiff_llhfit)
		VertexRes_linefit.Fill(vertexdiff_linefit)
		DirectionRes_linefit_vs_minr.Fill(minr,solidangle_linefit)
		DirectionRes_llhfit_vs_minr.Fill(minr,solidangle_llhfit)


fout.cd()
DirectionRes_linefit.Write()
DirectionRes_llhfit.Write()
VertexRes_llhfit.Write()
VertexRes_linefit.Write()
DirectionRes_linefit_vs_minr.Write()
DirectionRes_llhfit_vs_minr.Write()
fout.Close()










