import numpy as np
from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate
import tables as tb

pdf = [[]]
pdfnorm = []
time_lim = [0.0,0.0]
dist_lim = [0.0,0.0]
peaktime = [0.0]

def readTables() :
	global pdf
	global time_lim
	global dist_lim
	global peaktime

	infile = open("fittertables.dat","r")
	lines = infile.readlines()
	linecount = 0
	xcount = 0
	ny = 0
	nx = 0
	minx = 0.0
	maxx = 0.0
	miny = 0.0
	maxy = 0.0
	maxvalue = 0.0;

	for line in lines :
		splitline = line.split(",",100)
		if linecount == 0 :
			nx = int(splitline[0].replace("\n",""))
			ny = int(splitline[1].replace("\n",""))
			minx = float(splitline[2].replace("\n",""))
			maxx = float(splitline[3].replace("\n",""))
			miny = float(splitline[4].replace("\n",""))
			maxy = float(splitline[5].replace("\n",""))
			linecount += 1
		else :
			if xcount == ny :
				pdf.append([])
				maxvalue = 0.0
				peaktime.append(0.0)
				xcount = 0
			for value in splitline :
				pdf[-1].append(float(value.replace("\n","")))
				if pdf[-1][-1] > maxvalue :
					maxvalue = pdf[-1][-1]
					peaktime[-1] = len(pdf[-1])-1
				xcount += 1
			linecount += 1
	time_lim = [miny,maxy]
	dist_lim = [minx,maxx]
	print(len(pdf))

def ComputeTotalChargeIntegrals() :
	global pdf
	global pdfnorm

	for i in range(len(pdf)) :
		integral = 0.0
		for j in range(len(pdf[i])) :
			integral += pdf[i][j]
		pdfnorm.append(integral)

def ComputeChiSqr(n=1.0,L=50.0) :
	global pdfnorm

	chi2 = 0.0;
	for d in range(10,100) :
		mc = pdfnorm[d]
		expected = n*np.exp(-(float(d)/L))	
		chi2 += ((mc-expected)**2.0)/mc

	return chi2

if __name__ == "__main__":

	readTables()
	ComputeTotalChargeIntegrals()

	n = 2.6071109198170515
	L = 18.949132224466762

	fout = ROOT.TFile("AttenuationPlots.root","RECREATE")
	MCAtten = ROOT.TH1F("MCAtten","",200,0.,200.)
	FitAtten = ROOT.TH1F("FitAtten","",200,0.0,200.)

	for i in range(10,200) :
		MCAtten.Fill(float(i),pdfnorm[i])
		FitAtten.Fill(float(i),n*np.exp(-float(i)/L)/(float(i)**2.0))

	fout.cd()
	MCAtten.Write()
	FitAtten.Write()
	fout.Close()
