from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",  type=str, default=None, help="csv file to read")
parser.add_argument("-o","--output", type=str, default="out", help="i3 file name to write")
parser.add_argument("-d","--ndoms",  type=int, default=20, help="Doms per string.")
parser.add_argument("-r","--domradius", type=int, default=(17.0 * 2.54 * 0.01 * 0.5), help='Radius of dom. Defaults to 17"')
parser.add_argument("-p", "--npmts", type=int, default=16, help="PMTs per DOM.")
args = parser.parse_args()


outfileName = (
    "PONE_" + str(args.output)+".i3.gz"
)

outfile = dataio.I3File(outfileName, "w")
domsPerString = args.ndoms

def generateGeometry():
    '''
    Creates an I3OMGeoMap object to return, created from x,y coordinates in a csv file

            Parameters:
                    None
            Returns:
                    I3OMGeoMap: map of Optical Modules (OMs)
    '''
    global domsPerString
    global spacing

    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0) #not sure what this does, but it was in other files
    area = 4.0 * ((args.domradius) ** 2.0) * np.pi * I3Units.meter2 
    geomap = dataclasses.I3OMGeoMap()


    #create list of depths for modules, currently hardcoded. Should add to params
    sp = 950.0 / 19.0
    depthlist = [(-450.0 + sp * i) * I3Units.meter for i in range(20)]
    depth = np.array(depthlist) 


    #read x,y from csv file
    xpositions = []
    ypositions = []

    if(args.input==None): raise Exception("no csv input file specified") 


    with open(args.input, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header if there is one
        for row in reader:
            xpositions.append(float(row[0]))
            ypositions.append(float(row[1]))

    #Add OMs at x,y and depth
    for i in range(len(xpositions)):
        for m in range(domsPerString):
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(
                xpositions[i], ypositions[i], depth[m]
            )
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
