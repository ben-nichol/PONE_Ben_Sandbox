#constructing STRAW geometry

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
from icecube.dataclasses import I3Constants
import numpy as np
import argparse
import sys
#sys.path.insert(0,'/home/users/dhilu/P_ONE_dvirhilu/src')

import gcdHelpers

outfileName = "STRAW_DOM_at_zero.i3.gz"
outfile = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/' + outfileName, 'w')

def generateGeometry():
    orientation = dataclasses.I3Orientation(0, 0, 1, 1, 0, 0)
    area = 0.5857538*I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    #radius = np.arange(37, 0, -(37/nCircles))
    #stringSpacing = 200/DPS
    #depth = np.arange(200, 0, -(200/DPS)) * I3Units.meter
    depth = np.array([53.21])
    xPos = ([0])
    yPos = ([0])

    #for i in range(0, len(radius)):
        #spacing = (2*np.pi*radius[i])/stringsPerCircle[i]
        #thetaDiff = spacing/radius[i]
        #theta = np.arange(0, 2*np.pi, thetaDiff)
        #x = radius[i] * np.cos(theta) * I3Units.meter
        #y = radius[i] * np.sin(theta) * I3Units.meter
        #xPos = np.append(xPos,x)
        #yPos = np.append(yPos,y)

    #array = np.ones(xPos.shape)
    #zPos = [array*depth[j] for j in range(0, domsPerString)]
    zPos = ()

    omkey = OMKey(1, 1, 0)
    omGeometry = dataclasses.I3OMGeo()
    omGeometry.orientation = orientation
    omGeometry.area = area
    omGeometry.position = dataclasses.I3Position(0, 0, 0)
    geomap[omkey] = omGeometry

    return geomap

geometry = dataclasses.I3Geometry()

geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = generateGeometry()

gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

gframe["I3Geometry"] = geometry

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
