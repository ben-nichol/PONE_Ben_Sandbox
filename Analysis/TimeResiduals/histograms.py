#check histograms

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib, csv
from icecube.phys_services import I3Calculator

def timeRes(file):
    pfile = dataio.I3File(str(file))
    x, y, z, wlen, time = ([]), ([]), ([]), ([]), ([])
    frame = pfile.pop_frame()
    az = ([])
    zen = ([])
    scat = ([])
    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                if photon.wavelength / I3Units.nanometer >= 400 and photon.wavelength / I3Units.nanometer <= 410:
                    photonPos = photon.pos
                    dir = photon.startDir
                    x = np.append(x, photonPos.x)
                    y = np.append(y, photonPos.y)
                    z = np.append(z, photonPos.z)
                    distance = np.sqrt(photonPos.x**2 + photonPos.y**2 + (photonPos.z - 107.66)**2)
                    if distance > 51 and distance < 55:
                       az = np.append(az, dir.azimuth)
                       zen = np.append(zen, dir.zenith)
                       wlen = np.append(wlen, photon.wavelength)
                       time = np.append(time, photon.time)
                       scat = np.append(scat, photon.numScattered)

    print("length of data - ", len(x))
    #distance = np.sqrt(x**2 + y**2 + (z - 107.66)**2)
    #boolean = [(distance > 51) & (distance < 55)]
    selectDist = distance#[boolean]
    selectTime = time#[boolean]

    return distance, az, zen, scat

dist, az, zen, scatters = ([]), ([]), ([]), ([])

for i in range(0,10):
    infile = '/data/p-one/akatil/timeResiduals/genPhotons_' + str(i) + '_violet.i3.gz'
    #infile = '/home/users/akatil/P-ONE/medium/SCPhotonsSTRAW.i3'
    d, a, z, s = timeRes(infile)
    dist = np.append(dist, d)
    az = np.append(az, a)
    zen = np.append(zen, z)
    scatters = np.append(scatters, s)

bins = np.linspace(50,55,301)
plt.figure()
plt.hist(dist, bins = bins,  histtype = 'step')
plt.xlabel("Time(ns)")
plt.ylabel("Count")
plt.title("Distance Distribution")
plt.show()

plt.savefig('DistanceDistribution.pdf', dpi = 200)
plt.clf()

plt.figure()
plt.hist(az, bins = 100,  histtype = 'step', label = 'Simulation')
plt.xlabel("angle")
plt.ylabel("Count")
plt.title("Azimuthal Distribution")
plt.show()

plt.savefig('AzimuthalDistribution.pdf', dpi = 200)
plt.clf()

plt.figure()
plt.hist((zen), bins = 100,  histtype = 'step', log=True, label = 'Simulation')
plt.xlabel("angle")
plt.ylabel("Count")
plt.title("Zenith Distribution")
plt.show()

plt.savefig('ZenithDistribution.pdf', dpi = 200)
plt.clf()

plt.figure()
plt.hist(scatters, bins = 100,  histtype = 'step', label = 'Simulation')
plt.xlabel("number of scatters")
plt.ylabel("Count")
plt.title("Scatter Distribution")
plt.show()

plt.savefig('ScatterDistribution.pdf', dpi = 200)
plt.clf()
