#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
args = parser.parse_args()

# open file
infileList = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep2/HorizGeo/NuGen_step2_HorizGeo_' + str(i) + '.i3.gz')
    infileList.append(infile)

# get distribution of direction, position, and energy
zenith = []
azimuth = []
energy = []
xdir = []

for infile in infileList:
    while( infile.more() ):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        zenith.append(np.cos(primary.dir.zenith))
        azimuth.append(primary.dir.azimuth)
        energy.append(primary.energy)
        xdir.append(-primary.dir.x)

logE = np.log10(energy)
print len(logE)
plt.hist(logE, histtype = "step", log = True)
plt.title("Weighted Muon Log Energy Distribution")
plt.xlabel("LogE (Energy in GeV, log base 10)")

plt.figure()
plt.hist(azimuth, histtype = "step", log = True)
plt.title("Weighted Muon Angular Distribution (Azimuth)")
plt.xlabel("Cosine of the Azimuth Angle")

plt.figure()
plt.hist(zenith, histtype = "step", log = True)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(xdir, histtype = "step")
plt.title("Relative Angle to String Distribution")
plt.xlabel("Cosine of the Relative Angle (x component of direction)")

plt.figure()
plotOtuputs = plt.hist2d(logE,xdir, norm = matplotlib.colors.LogNorm())
plt.title("Energy and Relative Angle Distribution")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.colorbar(plotOtuputs[3])

plt.show()