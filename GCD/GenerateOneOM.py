from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
from icecube.dataclasses import I3Constants
import numpy as np
import argparse
import sys
#from GCD.GenerateLatticeStructure import generateLaticeSpots
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--domradius",
    type=int,
    default=(17.0 * 2.54 * 0.01 * 0.5),
    help='Radius of dom. Defaults to 17"',
)
parser.add_argument("-p", "--npmts", type=int, default=16, help="PMTs per DOM.")
args = parser.parse_args()


outfileName = 'one-om-gcd-origin.i3.gz'
outfile = dataio.I3File(outfileName, "w")
nstrings = 1
spacing = 10
domsPerString = 1


def generateGeometry():
    global Rows
    global domsPerString
    global spacing

    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4.0 * ((args.domradius) ** 2.0) * np.pi * I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    offset = np.pi * (1.0 / 6)
    anglediff = np.pi * (1.0 / 3)
    neighbourangles = [
        offset,
        anglediff + offset,
        2.0 * anglediff + offset,
        3.0 * anglediff + offset,
        4.0 * anglediff + offset,
        5.0 * anglediff + offset,
    ]

    stringposx = [0.0]
    stringposy = [0.0]
    #    stringposx = [0.0]
    #    stringposy = [0.0]
    #
    #    while len(stringposx) < nstrings :
    #
    #        minradius = 1000000.0
    #        minradstring = 0
    #        minradstringneighbours = 10000
    #        for i in range(len(stringposx)) :
    #                nneighbours = 0
    #                rad = np.sqrt((stringposx[i])**2.0+(stringposy[i])**2.0)
    #                for j in range(len(stringposx)):
    #                        if i==j :
    #                                continue
    #                        dist = np.sqrt((stringposx[j]-stringposx[i])**2.0+(stringposy[j]-stringposy[i])**2.0)
    #                        if dist<1.2 :
    #                                nneighbours += 1
    #                if nneighbours < len(neighbourangles) and rad <= minradius :
    #                        if rad < minradius :
    #                                minradius = rad;
    #                                minradstring = i
    #                                minradstringneighbours = nneighbours
    #                        elif nneighbours < minradstringneighbours :
    #                                minradius = rad
    #                                minradstring = i
    #                                minradstringneighbours = nneighbours
    #
    #        maxneighours = 0
    #        maxneighbourstring = 0
    #        for j in range(len(neighbourangles)) :
    #                newposx = stringposx[minradstring]+np.sin(neighbourangles[j])
    #                newposy = stringposy[minradstring]+np.cos(neighbourangles[j])
    #
    #                nneighbours = 0
    #                overlap = False
    #                for k in range(len(stringposx)) :
    #                        dist = np.sqrt((newposx-stringposx[k])**2.0+(newposy-stringposy[k])**2.0)
    #                        if dist < 0.8 :
    #                                nneighbours = 0
    #                                overlap = True
    #                        if dist<1.2 :
    #                                nneighbours += 1
    #                if nneighbours > maxneighours and not overlap:
    #                        maxneighours = nneighbours
    #                        maxneighbourstring = j
    #        stringposx.append(stringposx[minradstring]+np.sin(neighbourangles[maxneighbourstring]))
    #        stringposy.append(stringposy[minradstring]+np.cos(neighbourangles[maxneighbourstring]))
    #
    #    mean_x = sum(stringposx)/len(stringposx)
    #    mean_y = sum(stringposy)/len(stringposy)

    #sp = 950.0 / 19.0
    depthlist = [0. * I3Units.meter]

    #print(len(stringposx))

    #for i in range(len(stringposx)):
    #    print(str(stringposx[i] * spacing) + " , " + str(stringposy[i] * spacing))

    for i in range(len(stringposx)):
        for m in range(domsPerString):
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
            # x = stringposx[i]-mean_x
            # y = stringposy[i]-mean_y
            x = stringposx[i]
            y = stringposy[i]
            z = depthlist[m]
            omGeometry.position = dataclasses.I3Position(x * spacing, y * spacing, z)
            for j in range(args.npmts):
                omkey = OMKey(i + 1, m + 1, j + 1)
                geomap[omkey] = omGeometry

    return geomap


geometry = dataclasses.I3Geometry()

geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geomap = generateGeometry()
geometry.omgeo = geomap

gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

geomap = generateGeometry()

gframe["I3Geometry"] = geometry
gframe["I3OMGeoMap"] = geomap
modgeomap = dataclasses.I3ModuleGeoMap()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    module = dataclasses.I3ModuleGeo()
    module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
    module.orientation = geomap[dom].orientation
    module.pos = geomap[dom].position
    module.radius = np.sqrt(geomap[dom].area / (4.0 * np.pi))
    modgeomap[mkey] = module

gframe["I3ModuleGeoMap"] = modgeomap
subdetec = dataclasses.I3MapModuleKeyString()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    subdetec[mkey] = "Upgrade"

gframe["Subdetectors"] = subdetec

gframe["StartTime"] = gcdHelpers.start_time
gframe["EndTime"] = gcdHelpers.end_time

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
