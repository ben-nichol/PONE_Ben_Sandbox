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
    '-r',
    '--domradius',
    type=int,
    default=(17.0 * 2.54 * 0.01 * 0.5),
    help='Radius of dom. Defaults to 17"',
)
parser.add_argument('-p', '--npmts', type=int, default=16, help='PMTs per DOM.')
args = parser.parse_args()


outfileName = 'one-om-gcd-origin.i3.gz'
outfile = dataio.I3File(outfileName, 'w')
nstrings = 1
spacing = 10
domsPerString = 3


def generateGeometry():
    global Rows
    global domsPerString
    global spacing

    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4.0 * ((args.domradius) ** 2.0) * np.pi * I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()


    stringposx = [0.0]
    stringposy = [0.0]

    depthlist = [0. * I3Units.meter, 250. * I3Units.meter, -250. * I3Units.meter]


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

gframe['I3Geometry'] = geometry
gframe['I3OMGeoMap'] = geomap
modgeomap = dataclasses.I3ModuleGeoMap()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    module = dataclasses.I3ModuleGeo()
    module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
    module.orientation = geomap[dom].orientation
    module.pos = geomap[dom].position
    module.radius = np.sqrt(geomap[dom].area / (4.0 * np.pi))
    modgeomap[mkey] = module

gframe['I3ModuleGeoMap'] = modgeomap
subdetec = dataclasses.I3MapModuleKeyString()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    subdetec[mkey] = 'Upgrade'

gframe['Subdetectors'] = subdetec

gframe['StartTime'] = gcdHelpers.start_time
gframe['EndTime'] = gcdHelpers.end_time

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
