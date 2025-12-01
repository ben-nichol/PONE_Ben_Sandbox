from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
from icecube.dataclasses import I3Constants
import numpy as np
import argparse
import sys
from GCD.GenerateLatticeStructure import generateLatticeSpots
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument(
    "-s", "--spacing", type=float, default=80.0, help="Spacing for strings."
)
parser.add_argument("-d", "--ndoms", type=int, default=20, help="Doms per string.")
parser.add_argument(
    "-r",
    "--domradius",
    type=int,
    default=(17.0 * 2.54 * 0.01 * 0.5),
    help='Radius of dom. Defaults to 17"',
)
parser.add_argument("-p", "--npmts", type=int, default=16, help="PMTs per DOM.")
args = parser.parse_args()


outfileName = "PONE_Special6String_" + str(int(args.spacing)) + "Spacing.i3.gz"
outfile = dataio.I3File(outfileName, "w")
nstrings = 6
spacing = args.spacing
domsPerString = args.ndoms


def generateGeometry():
    global Rows
    global domsPerString
    global spacing

    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4.0 * ((args.domradius) ** 2.0) * np.pi * I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    offset = 0.0
    anglediff = np.pi * (2.0 / 5)
    neighbourangles = [
        offset,
        anglediff + offset,
        2.0 * anglediff + offset,
        3.0 * anglediff + offset,
        4.0 * anglediff + offset,
    ]
    stringposx = [
        0.0,
        np.cos(neighbourangles[0]),
        np.cos(neighbourangles[1]),
        np.cos(neighbourangles[2]),
        np.cos(neighbourangles[3]),
        np.cos(neighbourangles[4]),
    ]
    stringposy = [
        0.0,
        np.sin(neighbourangles[0]),
        np.sin(neighbourangles[1]),
        np.sin(neighbourangles[2]),
        np.sin(neighbourangles[3]),
        np.sin(neighbourangles[4]),
    ]

    theta = [0.0, 0.0, 0.0, 0.0, 0.0]

    sp = 950.0 / 19.0
    depthlist = [(-450.0 + sp * i) * I3Units.meter for i in range(20)]

    for i in range(len(stringposx)):
        for m in range(domsPerString):
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
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
