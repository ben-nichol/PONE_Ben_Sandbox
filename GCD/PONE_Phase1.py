from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
from icecube.dataclasses import I3Constants
import numpy as np
import argparse
import sys

import gcdHelpers

outfileName = "PONE_Phase1.i3.gz"
outfile = dataio.I3File(outfileName, 'w')
numberOfCircles = 2
domsPerString  = 20
stringsPerCircle = ([7, 3])

def generateGeometry(nCircles, DPS, strings):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 0.5857538*I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    radius = np.arange(200, 0, -(200/nCircles))
    #stringsPerCircle = radius*ratio
    stringSpacing = 800/DPS
    depth = np.arange((I3Constants.SurfaceElev - I3Constants.OriginElev - 1600), (I3Constants.SurfaceElev - I3Constants.OriginElev - 2400), -(800/DPS)) * I3Units.meter
    xPos = ([])
    yPos = ([])
    for i in range(0, len(radius)):
        spacing = (2*np.pi*radius[i])/strings[i]
        thetaDiff = spacing/radius[i]
        theta = np.arange(0, 2*np.pi, thetaDiff)
        x = radius[i] * np.cos(theta) * I3Units.meter
        y = radius[i] * np.sin(theta) * I3Units.meter
        xPos = np.append(xPos,x)
        yPos = np.append(yPos,y)

    array = np.ones(xPos.shape)
    zPos = [array*depth[j] for j in range(0, domsPerString)]

    for m in range(0, DPS):
        for n in range(0, len(xPos)):
            omkey = OMKey(n, m, 0)
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.orientation = orientation
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(xPos[n], yPos[n], zPos[m][n])
            geomap[omkey] = omGeometry

    return geomap

geometry = dataclasses.I3Geometry()

geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = generateGeometry(numberOfCircles, domsPerString, stringsPerCircle)

gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

gframe["I3Geometry"] = geometry

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
