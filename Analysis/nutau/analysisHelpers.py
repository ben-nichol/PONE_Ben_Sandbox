#python script

from icecube import dataclasses, dataio, simclasses, icetray
from icecube.icetray import OMKey
import matplotlib.pylab as plt
from mpl_toolkits import mplot3d
import numpy as np
from I3Tray import I3Units
from matplotlib.colors import LogNorm

def plotGeo(geoFile, plotTitle):
    file = dataio.I3File(str(geoFile))
    frame = file.pop_frame()
    geometry = frame["I3Geometry"]
    omgeo = geometry.omgeo

    x, y, z = [], [], []
    for i in omgeo.keys():
        oKey = omgeo.get(i)
        domPos = oKey.position
        x.append(domPos.x)
        y.append(domPos.y)
        z.append(domPos.z)

    X, Y, Z = np.array(x), np.array(y), np.array(z)
    #distance = np.sqrt(X**2 + Y**2 + (100-Z)**2)
    #print(distance)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, zdir='z', s=20, c='b',rasterized=True)
    zlabel = ax.set_zlabel('\n$z$ (m)',fontsize=16)
    ylabel = ax.set_ylabel('\n$y$ (m)',fontsize=16)
    xlabel = ax.set_xlabel('\n$x$ (m)',fontsize=16)
    plt.title(str(plotTitle), fontsize=18)

    plt.savefig(str(plotTitle)+'_geo.pdf', dpi = 200)
    plt.clf()
    file.close()

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=20, c='b',rasterized=True)
    ylabel = ax.set_ylabel('\n$y$ (m)',fontsize=16)
    xlabel = ax.set_xlabel('\n$x$ (m)',fontsize=16)
    plt.title(str(plotTitle), fontsize=18)

    plt.savefig(str(plotTitle)+'_geo_xy.pdf', dpi = 200)
    plt.clf()
    file.close()

def plotI3Photons(file, plotTitle):
    pfile = dataio.I3File(str(file))
    x, y, z = [], [], []
    frame = pfile.pop_frame()

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                photonPos = photon.pos
                x.append(photonPos.x)
                y.append(photonPos.y)
                z.append(photonPos.z)

    X, Y, Z = np.array(x), np.array(y), np.array(z)
    r = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(r, Z-100)

    print(len(x))
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, zdir='z', s=20, c='b',rasterized=True)
    zlabel = ax.set_zlabel('\n$z$ (m)',fontsize=16)
    ylabel = ax.set_ylabel('\n$y$ (m)',fontsize=16)
    xlabel = ax.set_xlabel('\n$x$ (m)',fontsize=16)
    plt.title(str(plotTitle), fontsize=18)
    plt.clf()

    np.savetxt(str(plotTitle) + '.csv', phi, delimiter=',')

    #plt.savefig(str(plotTitle)+'_hits.pdf', dpi = 200)
    plt.clf()
    pfile.close()

def histograms(file, selectDist, zenithDistribution):
    pfile = dataio.I3File(str(file))
    x, y, z, wlen, time = ([]), ([]), ([]), ([]), ([])
    frame = pfile.pop_frame()
    az = ([])
    zen = ([])
    scat = ([])
    startTime = ([])

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                #if photon.wavelength / I3Units.nanometer >= 400 and photon.wavelength / I3Units.nanometer <= 410:
                photonPos = photon.pos
                dir = photon.startDir
                x = np.append(x, photonPos.x)
                y = np.append(y, photonPos.y)
                z = np.append(z, photonPos.z)
                distance = np.sqrt(photonPos.x**2 + photonPos.y**2 + (photonPos.z - 107.66)**2)
                print (distance)
                #if distance < selectDist+3:  #distance > selectDist-3 and
                az = np.append(az, dir.azimuth)
                zen = np.append(zen, dir.zenith)
                wlen = np.append(wlen, photon.wavelength)
                time = np.append(time, photon.time)
                scat = np.append(scat, photon.numScattered)
                startTime = np.append(startTime, photon.startTime)
    plt.figure()
    plt.hist(az, bins = 100,  histtype = 'step', log=True, label = 'Simulation')
    plt.xlabel("angle")
    plt.ylabel("Count")
    plt.title("Azimuthal Distribution_" + str(selectDist) +'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('TestingFlasher_AzimuthalDistribution_' + str(int(selectDist)) + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    plt.figure()
    plt.hist(np.cos(zen), bins = 100,  histtype = 'step', log=True, label = 'Simulation')
    plt.xlabel("cos(angle)")
    plt.ylabel("Count")
    plt.title("cos(Zenith) Distribution_" + str(selectDist)+'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('TestingFlasher_CosZenithDistribution_' + str(int(selectDist)) + '_'+ str(zenithDistribution) + '.pdf', dpi = 200)
    plt.clf()

    plt.figure()
    plt.hist(scat, bins = 100,  histtype = 'step', label = 'Simulation')
    plt.xlabel("number of scatters")
    plt.ylabel("Count")
    plt.title("Scatter Distribution_" + str(selectDist) +'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('TestingFlasher_ScatterDistribution_' + str(int(selectDist)) + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    plt.figure()
    plt.hist2d(np.cos(zen), scat, bins = 100, norm = LogNorm())
    cb = plt.colorbar()
    #plt.plot(np.cos(zen), scat, '.')
    plt.xlabel("cos(zen)")
    plt.ylabel("NumScattered")
    plt.title("Scattering vs cos(Zenith)"+'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('ScattervsCosZenHist' + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    plt.figure()
    plt.hist2d(az, scat, bins = 100, norm = LogNorm())
    cb = plt.colorbar()
    #plt.plot(np.cos(zen), scat, '.')
    plt.xlabel("azimuth")
    plt.ylabel("NumScattered")
    plt.title("Scattering vs Azimuth"+'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('ScattervsAzimuthHist' + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    plt.figure()
    plt.hist(startTime, bins = 100,  histtype = 'step', label = 'Simulation')
    plt.xlabel("StartTime(ns)")
    plt.ylabel("Count")
    plt.title("StartTime Distribution_" + str(selectDist) +'_'+ str(zenithDistribution))
    plt.show()

    plt.savefig('TestingFlasher_StartTimeDistribution_' + str(int(selectDist)) + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    pfile.close()



def csvI3PhotonsPos(file, filename, zenithDistribution):
    pfile = dataio.I3File(str(file))
    x, y, z = ([]), ([]), ([])
    frame = pfile.pop_frame()

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                photonPos = photon.pos
                x = np.append(x, photonPos.x)
                y = np.append(y, photonPos.y)
                z = np.append(z, photonPos.z)

    combined = np.vstack((x, y, z)).T

    np.savetxt(str(filename) +'_'+ str(zenithDistribution) + '.csv', combined, delimiter=',')
    #np.savetxt(str(filename) +'_'+ str(zenithDistribution) + '.csv', combined)
    pfile.close()

def Dim2histograms(file, zenithDistribution):
    pfile = dataio.I3File(str(file))
    x, y, z, wlen, time = ([]), ([]), ([]), ([]), ([])
    frame = pfile.pop_frame()
    az = ([])
    zen = ([])
    scat = ([])
    startTime = ([])
    Zenith180PhotonCount = 0
    Zenith42PhotonCount = 0
    Zenith90PhotonCount = 0

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                #if photon.wavelength / I3Units.nanometer >= 400 and photon.wavelength / I3Units.nanometer <= 410:
                photonPos = photon.pos
                dir = photon.startDir
                x = np.append(x, photonPos.x)
                y = np.append(y, photonPos.y)
                z = np.append(z, photonPos.z)
                #distance = np.sqrt(photonPos.x**2 + photonPos.y**2 + (photonPos.z - 107.66)**2)
                #print (distance)
                #if distance < selectDist+3:  #distance > selectDist-3 and
                az = np.append(az, dir.azimuth)
                zen = np.append(zen, dir.zenith)
                wlen = np.append(wlen, photon.wavelength)
                time = np.append(time, photon.time)
                scat = np.append(scat, photon.numScattered)
                startTime = np.append(startTime, photon.startTime)

    plt.figure()
    plt.hist2d(np.cos(zen), az, bins = 1000)
    cb = plt.colorbar()
    plt.xlabel("cos(zen)")
    plt.ylabel("azi")
    plt.title("AzivsCosZen" + str(zenithDistribution))
    plt.show()

    plt.savefig('AzivsCosZen' + '_'+ str(zenithDistribution) +'.pdf', dpi = 200)
    plt.clf()

    return az, zen
