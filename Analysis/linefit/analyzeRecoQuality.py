#!/usr/bin/env python

# Modified to make plots that compare and analyze the reconstruction data. -dg

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import statistics as stats
from SimAnalysis import makeRatioHist
import argparse

parser = argparse.ArgumentParser(description = "Compares reconstrution of the muons to the simulated tracks")
parser.add_argument( '-D', '--domthresh', default = "7", help = "number of doms needed to be hit for a frame to be considered")
parser.add_argument( '-t', '--recotype', help = "type of reconstruction (linefit/improved)")
parser.add_argument( '-g', '--gcd', help = "gcd used")
parser.add_argument( '-r', '--runnum', default = "1", help = "Run number used to generate file with NuGen")
parser.add_argument( '-n', '--num', default = "10", help = "Number of Energy bins to use for plotting purposes")

args = parser.parse_args()

# Use the type of file to determine where it is
if args.recotype == "linefit":
    direc = "linefitReco"
elif args.recotype == "improved":
    direc = "improvedfitReco"
elif args.recotype == "llhfit":
    direc = "llhfitReco"
else:
    print("Unrecognized Type! Currently accepts 'improved', 'linefit' or 'llhfit'.")
    exit()

infile = dataio.I3File("/home/users/ghuman/simAnalysis/output/I3Files/" + direc + "/" + args.gcd + "/NuGen_" + direc + "_" + args.gcd + "_d" + str(args.domthresh) + "_14" + str(args.runnum) + ".i3.gz")

def findIndex(energy, bins):
    logE = np.log10(energy)

    for i in range(len(bins)):
        if logE < bins[i+1]:
            return i
    
    raise ValueError("energy not in range")

def fitanalysis(infile, binsE, fitType, n = 10):
    if fitType == "llhfit":
        cosAlpha = [[],[],[]]
        listsError = [[[] for i in range(n)],[[] for i in range(n)],[[] for i in range(n)]]
    elif fitType == "improved":
        cosAlpha = [[],[]]
        listsError = [[[] for i in range(n)],[[] for i in range(n)]]
    else:
        cosAlpha = []
        listsError = [[] for i in range(n)]
    unsuccessfulRecoEnergy = []
    successfulRecoEnergy = []
    energy = []
    orgRelSpeed = [[] for i in range(n)]
    relSpeed = []

    for frame in infile:
        primary = frame["NuGPrimary"]
        mctree = frame["I3MCTree"]
        muon = dataclasses.I3MCTree.first_child(mctree, primary)
        
        if fitType == "linefit":
            recoParticle = frame["LineFitRecoParticle"]
        elif fitType == "improved":
            recoParticle = frame["ImprovedRecoParticle"]
            recoParticle2 = frame["LineFitRecoParticle"]
        else:
            recoParticle = frame["LlhFitRecoParticle"]
            recoParticle2 = frame["LineFitRecoParticle"]

        if recoParticle.fit_status == dataclasses.I3Particle.InsufficientQuality:
            unsuccessfulRecoEnergy.append(np.log10(primary.energy))
            continue

        muonDir = muon.dir
        recoDir = recoParticle.dir
        index = findIndex(primary.energy, binsE)

        dotProduct = muonDir.x*recoDir.x + muonDir.y*recoDir.y + muonDir.z*recoDir.z

        if fitType == "linefit":
            relSpeed.append(abs(recoParticle.speed - muon.speed)/muon.speed)
            orgRelSpeed[index].append(abs(recoParticle.speed - muon.speed)/muon.speed)
        else:
            relSpeed.append(abs(recoParticle2.speed - muon.speed)/muon.speed)
            orgRelSpeed[index].append(abs(recoParticle2.speed - muon.speed)/muon.speed)

        if fitType != "linefit":
            recoDir2 = recoParticle2.dir
            dotProduct2 = muonDir.x*recoDir2.x + muonDir.y*recoDir2.y + muonDir.z*recoDir2.z
            cosAlpha[0].append(dotProduct)
            cosAlpha[1].append(dotProduct2)
            listsError[0][index].append(np.arccos(dotProduct)/I3Units.deg)
            listsError[1][index].append(np.arccos(dotProduct2)/I3Units.deg)
        else:
            cosAlpha.append(dotProduct)
            listsError[index].append(np.arccos(dotProduct)/I3Units.deg)   

        successfulRecoEnergy.append(np.log10(primary.energy))
        energy.append(np.log10(primary.energy))

    meanRelSpeed = []
    errorRelSpeed = []

    for speed in orgRelSpeed:
        if len(speed) < 2:            
            errorRelSpeed.append(0)
        else:
            errorRelSpeed.append(stats.stdev(speed))

        if len(speed) < 1:
            meanRelSpeed.append(0)
        else:
            meanRelSpeed.append(stats.mean(speed))

    if fitType == "linefit":
        alpha = [np.arccos(cosA)/I3Units.deg for cosA in cosAlpha]
        percent50Error = []
        percent90Error = []
        stddev = []
        mean = []
        for errorList in listsError:
            errorList.sort()
            index50Per = int(0.5*len(errorList))
            index90Per = int(0.9*len(errorList))
    
            if len(errorList) == 0:
                percent50Error.append(-10)
                percent90Error.append(-10)
                mean.append(-10)
                stddev.append(0)
            elif len(errorList) == 1:
                percent50Error.append(errorList[index50Per])
                percent90Error.append(errorList[index90Per])
                stddev.append(0)
                mean.append(stats.mean(errorList))
            else:
                percent50Error.append(errorList[index50Per])
                percent90Error.append(errorList[index90Per])
                stddev.append(stats.stdev(errorList))
                mean.append(stats.mean(errorList))

    else :
        alpha = [[np.arccos(cosA)/I3Units.deg for cosA in cosAlpha[0]],[np.arccos(cosA)/I3Units.deg for cosA in cosAlpha[1]]]
        percent50Error = [[],[]]
        percent90Error = [[],[]]
        stddev = [[],[]]
        mean = [[],[]]
        for i in range(len(listsError[0])):
            errorList = listsError[0][i]
            errorList2 = listsError[1][i]
            errorList.sort()
            errorList2.sort()
            index50Per = int(0.5*len(errorList))
            index50Per2 = int(0.5*len(errorList2))
            index90Per = int(0.9*len(errorList))
            index90Per2 = int(0.9*len(errorList2))
    
            if len(errorList) == 0:
                percent50Error[0].append(-10)
                percent90Error[0].append(-10)
                mean[0].append(-10)
                stddev[0].append(0)

            elif len(errorList) == 1:
                percent50Error[0].append(errorList[index50Per])
                percent90Error[0].append(errorList[index90Per])
                stddev[0].append(0)
                mean[0].append(stats.mean(errorList))
            else:
                percent50Error[0].append(errorList[index50Per])
                percent90Error[0].append(errorList[index90Per])
                stddev[0].append(stats.stdev(errorList))
                mean[0].append(stats.mean(errorList))

            if len(errorList2) == 0:
                percent50Error[1].append(-10)
                percent90Error[1].append(-10)
                mean[1].append(-10)
                stddev[1].append(0)

            elif len(errorList2) == 1:
                percent50Error[1].append(errorList2[index50Per2])
                percent90Error[1].append(errorList2[index90Per2])
                stddev[1].append(0)
                mean[1].append(stats.mean(errorList2))
            else:
                percent50Error[1].append(errorList2[index50Per2])
                percent90Error[1].append(errorList2[index90Per2])
                stddev[1].append(stats.stdev(errorList2))
                mean[1].append(stats.mean(errorList2))

    return percent50Error, percent90Error, alpha, cosAlpha, unsuccessfulRecoEnergy, successfulRecoEnergy, mean, stddev, energy, [entry*100 for entry in relSpeed], meanRelSpeed, errorRelSpeed


def main():
    
    binsE = np.linspace(2,10,int(args.num) + 1)
    fitAnalysis = fitanalysis(infile, binsE, args.recotype, int(args.num))

    # energy distribution works independent of fit type
    plt.figure()
    plt.hist(fitAnalysis[8], log = True, histtype = 'step', bins = 20, label = args.recotype)
    plt.xlabel(r'$log_{10}\, E/GeV$')
    plt.title('Energy Distribution in simulated events')
    plt.savefig('/home/users/ghuman/simAnalysis/output/plots/'+ direc + '/pentagon10040/energy_dist.png',dpi=300)

    # Absolute relative Speed distribution is purely from linefit, since improved fixes speed
    plt.figure()
    plt.hist(fitAnalysis[9], log = True, histtype = 'step', bins = 20, label = args.recotype)
    plt.xlabel(r'Percent Difference in Speed')
    plt.title('Percent Speed Difference Distribution')
    plt.savefig('/home/users/ghuman/simAnalysis/output/plots/'+ direc + '/pentagon10040/speed_dist.png',dpi=300)

    # Energy vs Mean angular error with stdev as errorbar
    plt.figure()
    plt.errorbar(binsE[:-1], fitAnalysis[10], yerr=fitAnalysis[11], fmt='.k', label = args.recotype)
    plt.xlabel(r'$log_{10}\, E/GeV$')
    plt.ylabel(r'Mean $(v_{i} - v)/v$')
    plt.title("Mean Relative Speed Difference")
    plt.legend()
    plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/meanSpeed_vs_energy.png',dpi=300)

    if args.recotype != "linefit":
        # Same plots as the simple linefit, but now with overlap

        # angle distribution
        plt.figure()
        plt.hist(fitAnalysis[2][0], log = True, histtype = 'step', bins = 20, label = args.recotype)
        plt.hist(fitAnalysis[2][1], log = True, histtype = 'step', bins = 20, label = "linefit")
        plt.xlabel(r'$\alpha$')
        plt.title('Distribution of Relative Angle of Muon and its Reconstruction')
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/'+ direc + '/pentagon10040/angle_dist.png',dpi=300)

        # Energy vs Mean angular error with stdev as errorbar
        plt.figure()
        eb1 = plt.errorbar(binsE[:-1], fitAnalysis[6][0], yerr=fitAnalysis[7][0], fmt='.k', label = args.recotype)
        eb2 = plt.errorbar(binsE[:-1], fitAnalysis[6][1], yerr=fitAnalysis[7][1], fmt='.g', label = "linefit")
        eb1[-1][0].set_linestyle('--')
        eb2[-1][0].set_linestyle('--')
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Error (degrees)")
        plt.title("Average Reconstruction Error")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/meanerror_vs_energy.png',dpi=300)

        # 50% of recos in this energy have a smaller error
        plt.figure()
        plt.step(binsE[:-1], fitAnalysis[0][0], where = 'post', label = args.recotype)
        plt.step(binsE[:-1], fitAnalysis[0][1], where = 'post', label = "linefit")
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Difference (degrees)")
        plt.title("Reconstuction Error - Successful Reco Only, 50th Percentile")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/50recoerror.png',dpi=300)
    
        # 90% of recos in this energy have a smaller error
        plt.figure()
        plt.step(binsE[:-1], fitAnalysis[1][0], where = 'post', label = args.recotype)
        plt.step(binsE[:-1], fitAnalysis[1][1], where = 'post', label = "linefit")
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Difference (degrees)")
        plt.title("Reconstuction Error - Successful Reco Only, 90th Percentile")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/90recoerror.png',dpi=300)

    else:
        # angle distribution
        plt.figure()
        plt.hist(fitAnalysis[2], log = True, histtype = 'step', bins = 20, label = args.recotype)
        plt.xlabel(r'$\alpha$')
        plt.title('Distribution of Relative Angle of Muon and its Reconstruction')
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/'+ direc + '/pentagon10040/angle_dist.png',dpi=300)

        # Energy vs Mean angular error with stdev as errorbar
        plt.figure()
        plt.errorbar(binsE[:-1], fitAnalysis[6], yerr=fitAnalysis[7], fmt='.k', label = args.recotype)
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Error (degrees)")
        plt.title("Average Reconstruction Error")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/meanerror_vs_energy.png',dpi=300)

        # 50% of recos in this energy have a smaller error
        plt.figure()
        plt.step(binsE[:-1], fitAnalysis[0], where = 'post', label = args.recotype)
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Difference (degrees)")
        plt.title("Reconstuction Error - Successful Reco Only, 50th Percentile")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/50recoerror.png',dpi=300)
    
        # 90% of recos in this energy have a smaller error
        plt.figure()
        plt.step(binsE[:-1], fitAnalysis[1], where = 'post', label = args.recotype)
        plt.xlabel(r'$log_{10}\, E/GeV$')
        plt.ylabel("Angular Difference (degrees)")
        plt.title("Reconstuction Error - Successful Reco Only, 90th Percentile")
        plt.legend()
        plt.savefig('/home/users/ghuman/simAnalysis/output/plots/' + direc + '/pentagon10040/90recoerror.png',dpi=300)

'''
binsE = np.linspace(3,7,10)

plt.figure()
ratioHist, edges = makeRatioHist(fitanalysisImproved[4], fitanalysisImproved[5], weights1 = np.ones(len(fitanalysisImproved[4])), weights2 = np.ones(len(fitanalysisImproved[5])), bins = binsE)
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Fraction")
plt.title("Fraction of Failed Reconstructions")
'''
#plt.show()

if __name__ == "__main__":
    main()
