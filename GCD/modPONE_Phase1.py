#!/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/bin/python
from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
from icecube.dataclasses import I3Constants
import numpy as np
import argparse
import sys

# sys.path.insert(0,'/home/users/dhilu/P_ONE_dvirhilu/src')

import gcdHelpers

outfileName = "SingleOM17.i3.gz"
outfile = dataio.I3File(outfileName, "w")
numberOfCircles = 1
domsPerString = 1
stringsPerCircle = [3]


def generateGeometry(nCircles, DPS, strings):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4 * np.pi * (0.2159**2) * I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    radius = np.arange(500, 0, -(500 / nCircles))
    # stringsPerCircle = radius*ratio
    stringSpacing = 800 / DPS
    depth = (
        np.array([0]) * I3Units.meter
    )  # np.arange((I3Constants.SurfaceElev - I3Constants.OriginElev - 1600), (I3Constants.SurfaceElev - I3Constants.OriginElev - 2400), -(800/DPS)) * I3Units.meter
    xPos = []
    yPos = []
    for i in range(0, len(radius)):
        print("Making x and y positions")
        spacing = (2 * np.pi * radius[i]) / strings[i]
        thetaDiff = spacing / radius[i]
        theta = np.arange(0, 2 * np.pi, thetaDiff)
        x = radius[i] * np.cos(theta) * I3Units.meter
        y = radius[i] * np.sin(theta) * I3Units.meter
        xPos = np.append(xPos, x)
        yPos = np.append(yPos, y)

    array = np.ones(xPos.shape)
    zPos = [array * depth[j] for j in range(0, domsPerString)]

    for m in range(0, DPS):
        for n in range(0, len(xPos)):
            print("Placeing OMs")
            print(n, m, 0)
            omkey = OMKey(n, m, 0)
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.orientation = orientation
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(xPos[n], yPos[n], zPos[m][n])
            geomap[omkey] = omGeometry

    print("Placing Final OM")
    omkey = OMKey(3, 0, 0)
    omGeometry = dataclasses.I3OMGeo()
    omGeometry.orientation = orientation
    omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
    omGeometry.area = area
    omGeometry.position = dataclasses.I3Position(0, 0, 0)
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
gframe["I3OMGeoMap"] = geometry.omgeo
modgeomap = dataclasses.I3ModuleGeoMap()
for dom in geometry.omgeo.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    module = dataclasses.I3ModuleGeo()
    module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
    module.orientation = geometry.omgeo[dom].orientation
    module.pos = geometry.omgeo[dom].position
    module.radius = np.sqrt(geometry.omgeo[dom].area / (4.0 * np.pi))
    modgeomap[mkey] = module

gframe["I3ModuleGeoMap"] = modgeomap
subdetec = dataclasses.I3MapModuleKeyString()
for dom in geometry.omgeo.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    subdetec[mkey] = "Upgrade"

gframe["Subdetectors"] = subdetec

gframe["StartTime"] = gcdHelpers.start_time
gframe["EndTime"] = gcdHelpers.end_time

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
