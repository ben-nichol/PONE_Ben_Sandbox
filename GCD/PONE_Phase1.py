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
NPMTs = 13

def generateGeometry(nCircles, DPS, strings):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4.0*((17.0*2.54*0.01*0.5)**2.0)*np.pi*I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    radius = np.arange(200, 0, -(200/nCircles))
    #stringsPerCircle = radius*ratio
    #stringSpacing = 800/DPS
    #depth = np.arange(-500,500, -(1000/DPS)) * I3Units.meter
    sp = 800.0/19.0
    depthlist = [(-400.0+sp*i)*I3Units.meter for i in range(20)]
    depth = np.array(depthlist)
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
            #omkey = OMKey(n+1, m+1, 0)
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(xPos[n], yPos[n], zPos[m][n])
            for i in range(NPMTs) :
                omkey = OMKey(n+1, m+1, i+1)
                geomap[omkey] = omGeometry

    return geomap

geometry = dataclasses.I3Geometry()

geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geomap = generateGeometry(numberOfCircles, domsPerString, stringsPerCircle)
geometry.omgeo = geomap

gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

geomap = generateGeometry(numberOfCircles, domsPerString, stringsPerCircle)

gframe["I3Geometry"] = geometry
gframe["I3OMGeoMap"] = geomap
modgeomap = dataclasses.I3ModuleGeoMap()
for dom in geomap.keys() :
        mkey = dataclasses.ModuleKey(dom.string,dom.om)
        module = dataclasses.I3ModuleGeo()
        module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
        module.orientation = geomap[dom].orientation
        module.pos = geomap[dom].position
        module.radius = np.sqrt(geomap[dom].area/(4.0*np.pi))
        modgeomap[mkey] = module

gframe["I3ModuleGeoMap"] = modgeomap;
subdetec = dataclasses.I3MapModuleKeyString() 
for dom in geomap.keys() :
    mkey = dataclasses.ModuleKey(dom.string,dom.om)
    subdetec[mkey] = "Upgrade"
        
gframe["Subdetectors"] = subdetec

gframe["StartTime"] = gcdHelpers.start_time
gframe["EndTime"] = gcdHelpers.end_time

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
