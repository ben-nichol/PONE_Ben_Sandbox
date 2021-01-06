#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib
from simAnalysis import SimAnalysis

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
parser.add_argument( '-H', '--hitThresh', default = 1, help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-D', '--domThresh', default = 8, help = "threshold of hit DOMs for the frame to be considered")
parser.add_argument( '-R', '--maxResidual', default = 100 , help = "maximum time residual allowed for the hit to be considered")
args = parser.parse_args()

hitThresh = int(args.hitThresh)
domThresh = int(args.domThresh)
maxResidual = float(args.maxResidual)

#infileListOrig = []
#for i in range(300, 400):
#    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
#    infileListOrig.append(infile)

#gcdOrig = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz')
#geoMapOrig = gcdOrig.pop_frame()["I3Geometry"].omgeo

infileListNew = []
for i in range(500, 700):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/denseGeo/NuGen_step3_denseGeo_' + str(i) + '.i3.gz')
    infileListNew.append(infile)

#gcdl0 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_SL0.i3.gz')
gcdl1 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_SL1.i3.gz')
#gcdl2 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_SL2.i3.gz')
#geoMapl0 = gcdl0.pop_frame()["I3Geometry"].omgeo
geoMapl1 = gcdl1.pop_frame()["I3Geometry"].omgeo
#geoMapl2 = gcdl2.pop_frame()["I3Geometry"].omgeo

gcdTest = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_10LineGeometry.i3.gz') 
geoMapTest = gcdTest.pop_frame()["I3Geometry"].omgeo

#dataOrig, weightsOrig, _binsOrig = SimAnalysis.getEffectiveAreaData(infileListOrig, geoMapOrig.keys(), hitThresh, domThresh, maxResidual, geoMapOrig, 10)
#print "finished orig"
#datal0, weightsl0, _binsl1 = SimAnalysis.getEffectiveAreaData(infileListNew, geoMapl0.keys(), hitThresh, domThresh, maxResidual, geoMapl0, 10)
#print "finished l0"
#for infile in infileListNew:
#    infile.rewind()
datal1, weightsl1, binsl1 = SimAnalysis.getEffectiveAreaData(infileListNew, geoMapl1.keys(), hitThresh, domThresh, maxResidual, geoMapl1, 10)
print "finished l1"
#for infile in infileListNew:
#    infile.rewind()
#datal2, weightsl2, _binsl2 = SimAnalysis.getEffectiveAreaData(infileListNew, geoMapl2.keys(), hitThresh, domThresh, maxResidual, geoMapl2, 10)
#print "finished l2"
for infile in infileListNew:
    infile.rewind()
dataTest, weightsTest, _binsTest = SimAnalysis.getEffectiveAreaData(infileListNew, geoMapTest.keys(), hitThresh, domThresh, maxResidual, geoMapTest, 10)

binsE = binsl1["logEnergy"]
binsZenith = binsl1["cosZenith"]


plt.rcParams["mathtext.fontset"] = "cm"
'''
plt.figure()
plotOutputs = plt.hist2d(dataOrig["logEnergy"], dataOrig["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsOrig["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') Original')
plt.colorbar(plotOutputs[3])

plt.figure()
plotOutputs = plt.hist2d(datal0["logEnergy"], datal0["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl0["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') l0')
plt.colorbar(plotOutputs[3])
'''
plt.figure()
plotOutputs = plt.hist2d(datal1["logEnergy"], datal1["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl1["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') 5 Line Geometry')
plt.colorbar(plotOutputs[3])
'''
plt.figure()
plotOutputs = plt.hist2d(datal2["logEnergy"], datal2["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl2["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') l2')
plt.colorbar(plotOutputs[3])
'''
plt.figure()
plotOutputs = plt.hist2d(dataTest["logEnergy"], dataTest["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsTest["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') 10 Line Geometry')
plt.colorbar(plotOutputs[3])

plt.figure()
#plt.hist(dataOrig["logEnergy"], log = True, histtype = 'step', weights = weightsOrig["logEnergy"], bins = binsE, label = 'Orig')
#plt.hist(datal0["logEnergy"], histtype = 'step', log = True, weights = weightsl0["logEnergy"], bins = binsE, label = 'l0')
plt.hist(datal1["logEnergy"], histtype = 'step', log = True, weights = weightsl1["logEnergy"], bins = binsE, label = '5 Line Geometry', color = 'blue')
#plt.hist(datal2["logEnergy"], histtype = 'step', log = True, weights = weightsl2["logEnergy"], bins = binsE, label = '12')
plt.hist(dataTest["logEnergy"], histtype = 'step', log = True, weights = weightsTest["logEnergy"], bins = binsE, label = '10 Line geometry', color = 'violet')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("counts")
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.legend(loc = 'upper left')

plt.figure()
#plt.hist(dataOrig["Horizontal"], log = True, histtype = 'step', weights = weightsOrig["Horizontal"], bins = binsE, label = 'Orig')
#plt.hist(datal0["Horizontal"], histtype = 'step', log = True, weights = weightsl0["Horizontal"], bins = binsE, label = 'l0')
plt.hist(datal1["Horizontal"], histtype = 'step', log = True, weights = weightsl1["Horizontal"], bins = binsE, label = '5 Line Geometry', color = 'darkturquoise')
#plt.hist(datal2["Horizontal"], histtype = 'step', log = True, weights = weightsl2["Horizontal"], bins = binsE, label = '12')
plt.hist(dataTest["Horizontal"], histtype = 'step', log = True, weights = weightsTest["Horizontal"], bins = binsE, label = '10 Line geometry', color = 'violet')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("counts")
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') ' + r'$cos{\theta} \epsilon(-0.1,0.1)$')
plt.legend(loc = 'upper left')
'''
plt.figure()
ratioHist, edges = SimAnalysis.makeRatioHist(datal1["logEnergy"], dataTest["logEnergy"], weightsl1["logEnergy"], weightsTest["logEnergy"], binsE)
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of initial/new (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM)")

plt.figure()
ratioHist, edges = SimAnalysis.makeRatioHist(datal1["Horizontal"], dataTest["Horizontal"], weightsl1["Horizontal"], weightsTest["Horizontal"], binsE)
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of initial/new (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM), " + r'$cos{\theta} \epsilon(-0.1,0.1)$')
'''
plt.show()
