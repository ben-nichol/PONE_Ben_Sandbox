from icecube import dataclasses, dataio, simclasses, icetray
from icecube.phys_services import I3Calculator
from icecube.icetray import OMKey, I3Units, I3Frame
import matplotlib.pylab as plt
from mpl_toolkits import mplot3d
import numpy as np
from I3Tray import I3Units


gcd = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
geometry = gcd.pop_frame()["I3Geometry"]

directory='/data/p-one/akatil/step_2_medium_water/'
file = 'step_2_8_medium_water.i3.gz'

pfile = dataio.I3File(directory+str(file))
frame = pfile.pop_frame()

for pframe in pfile:
    x, y, z = ([]), ([]), ([])
    distance = ([])
    closestDistance = 1e7
    photonMap = pframe['I3Photons']
    primary = pframe["NuGPrimary"]
    mctree = pframe["I3MCTree"]
    tau = dataclasses.I3MCTree.first_child(mctree, primary)
    tauPosition = tau.pos
    for modKey,photons in photonMap:
        omkey = OMKey(modKey.string, modKey.om, 0)
        position = geometry.omgeo[omkey].position
        distance = np.append(distance, np.sqrt((tauPosition.x - position.x)**2 + (tauPosition.y - position.y)**2 + (tauPosition.z - position.z)**2))
        #tau.shape = dataclasses.I3Particle.Cascade
        #distance = I3Calculator.cherenkov_distance(tau, position)
        #if distance < closestDistance:
            #closestDistance = distance
            #closestOMPosition = position
        for photon in photons:
            photonPos = photon.pos
            x = np.append(x, photonPos.x)
            y = np.append(y, photonPos.y)
            z = np.append(z, photonPos.z)

    distanceBetween = np.sqrt((tauPosition.x - x)**2 + (tauPosition.y - y)**2 + (tauPosition.z - z)**2)
    difference = abs(distanceBetween-distance.min())
    print(tau.energy, distance.min(), len(distanceBetween[(difference < difference.min()+5) & (difference > difference.min()-5)]))
