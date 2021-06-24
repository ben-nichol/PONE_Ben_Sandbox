#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np

# open file
infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step1/MuonGun_step1_139005_000000.i3.bz2')
#outfile = open('muonDirectionz.txt', 'w+')

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())

# get distribution of direction and energy
zenith = []
azimuth = []
energy = []
zdirection = []
x = []
y = []
z = []
weight = []
#i=1
for frame in qframes:
    # get propagated muon
    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    muon = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

    # add muon energy
    energy.append(muon.energy / I3Units.GeV)
    weighter = frame["I3MCWeightDict"]
    weight.append(weighter["weight"])

    # get muon direction
    muonDirection = muon.dir
    azimuth.append(np.cos(muonDirection.azimuth))
    zenith.append(np.cos(muonDirection.zenith))
    #outfile.write("direction for muon number %d: \n" % (i))
    #outfile.write("x: " + str(muonDirection.x) + "\n")
    #outfile.write("y: " + str(muonDirection.y) + "\n")
    #outfile.write("z: " + str(muonDirection.z) + "\n")
    #outfile.write("azimuth: " + str(muonDirection.azimuth / I3Units.deg) + "\n")
    #outfile.write("zenith: " + str(muonDirection.zenith / I3Units.deg) + "\n")
    #outfile.write("phi: " + str(muonDirection.phi / I3Units.deg) + "\n")
    #outfile.write("theta: " + str(muonDirection.theta / I3Units.deg) + "\n\n\n")
    #i += 1

    # get muon position
    position = muon.pos
    x.append(position.x)
    y.append(position.y)
    z.append(position.z)


# plot histograms
xmin = 500 * I3Units.GeV
xmax = 1 * I3Units.TeV
plt.hist(energy, range = (xmin,xmax), histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Energy Distribution")
plt.xlabel("E(GeV)")
#plt.show()

plt.figure()
xmin = -1
xmax = 1
plt.hist(azimuth, range = (xmin,xmax), histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Azimuth)")
plt.xlabel("cos(phi)")
#plt.show()

plt.figure()
xmin = -1
xmax = 1
plt.hist(zenith, range = (xmin,xmax), histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel("Cosine of Zenith Angle")
#plt.show()

plt.figure()
plt.hist2d(zenith, energy, weights = weight, bins = 20)
plt.title("Weighted Distribution of Muons at Various Energies and Directions")
plt.xlabel("Cosine of Zenith Angle")
plt.ylabel("Energy (GeV)")
#plt.show()

plt.figure()
plt.hist(x, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (x)")
plt.xlabel("x Coordinate")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(y, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (y)")
plt.xlabel("y Coordinate")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(z, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (z)")
plt.xlabel("z Coordinate")
plt.ylabel("Number of Occurences")

plt.show()

#outfile.close()