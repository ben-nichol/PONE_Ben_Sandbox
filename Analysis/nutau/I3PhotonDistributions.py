from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib, csv
from statistics import mean
from icecube.phys_services import I3Calculator
import os

directory='/data/p-one/akatil/step_2_medium_water/'
file = ([])
for filename in os.listdir(directory):
    if filename.startswith("step_2"):
        if os.path.getsize(directory+str(filename)) > 0 and os.path.getsize(directory+str(filename)) < 2000000:
            file.append(filename)

x, y, z, az, zen, weights, wlen = ([]), ([]), ([]), ([]), ([]), ([]), ([])

#file = ['step_2_8_medium_water.i3.gz']
print(len(file))
for i in range(0, len(file)):
    print(i)
    pfile = dataio.I3File(directory+str(file[i]))
    frame = pfile.pop_frame()

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                photonPos = photon.pos
                photonDir = photon.dir
                x = np.append(x, photonPos.x)
                y = np.append(y, photonPos.y)
                z = np.append(z, photonPos.z)
                az = np.append(az, photonDir.azimuth)
                zen = np.append(zen, photonDir.zenith)
                weights = np.append(weights, photon.weight)
                wlen = np.append(wlen, photon.wavelength)


wScale = weights/weights.mean()
print("Weights Scaled - ", wScale.shape)
print("X - ", x.shape)

bins=1000
plt.hist(x, bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(x, bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend()
plt.xlabel("x(m)")
plt.title('I3Photon Position along x-Axis')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_xPositions.pdf', dpi=200)
plt.clf()

plt.hist(y, bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(y, bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend()
plt.xlabel("y(m)")
plt.title('I3Photon Position along y-Axis')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_yPositions.pdf', dpi=200)
plt.clf()

plt.hist(z, bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(z, bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend()
plt.xlabel("z(m)")
plt.title('I3Photon Position along z-Axis')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_zPositions.pdf', dpi=200)
plt.clf()

plt.hist(az, bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(az, bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend()
plt.xlabel("Azimuth")
plt.title('I3Photon Direction(Azimuth)')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_azimuth.pdf', dpi=200)
plt.clf()

plt.hist(np.cos(zen), bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(np.cos(zen), bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend()
plt.xlabel("Cos(Zenith)")
plt.title('I3Photon Direction(Zenith)')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_zenith.pdf', dpi=200)
plt.clf()

plt.hist(wlen, bins = bins, density = False, log=True, histtype='step', label = 'NoWeights')
plt.hist(wlen, bins = bins, density=False, log=True, histtype='step', weights=wScale, label='WithWeights')
plt.legend(loc='best')
plt.xlabel("Wavelength(m)")
plt.title('I3Photon Wavelength')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3Photons_wavelength.pdf', dpi=200)
plt.clf()
