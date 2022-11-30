from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
from GenerateLatticeStructure import generateLaticeSpots
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--spacing",type = float, default = 80.0, help="Spacing for strings in cluser.")
parser.add_argument("-l", "--clusterspacing",type = float, default = 400.0, help="Spacing for strings in cluser.")
parser.add_argument("-n", "--nstring",type= int, default = 10, help="Number of strings per cluster.")
parser.add_argument("-c", "--nclusters",type= int, default = 7, help="Number of clusters.")
parser.add_argument("-d", "--ndoms",type= int, default = 20, help="Doms per string.")
parser.add_argument("-r", "--domradius",type= int, default = (17.0*2.54*0.01*0.5), help="Radius of dom. Defaults to 17\"")
parser.add_argument("-p", "--npmts",type= int, default = 16, help="PMTs per DOM.")
args = parser.parse_args()


outfileName = "PONE_"+str(args.nstring)+"String_"+str(args.nclusters)+"Cluster.i3.gz"
outfile = dataio.I3File(outfileName, 'w')
nstrings = args.nstring
spacing = args.spacing
clusterspacing = args.clusterspacing
nclusters = args.nclusters
domsPerString  = args.ndoms

def generateGeometry():
    global Rows
    global domsPerString
    global spacing
    global clusterspacing

    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    area = 4.0*((args.domradius)**2.0)*np.pi*I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()


    stringposx, stringposy, theta = generateLaticeSpots(nstrings)

    clusterposx, clusterposy, clustertheta = generateLaticeSpots(nclusters)

    FinalStringx = []
    FinalStringy = []

    for i in range(len(clusterposx)) :
        for j in range(len(stringposx)):
            if i == 0 : 
                FinalStringx.append(stringposx[j]*spacing + clusterposx[i]*clusterspacing)
                FinalStringy.append(stringposy[j]*spacing + clusterposy[i]*clusterspacing)
            else :
                FinalStringx.append(np.cos(clustertheta[i])*stringposx[j]*spacing - np.sin(clustertheta[i])*stringposy[j]*spacing + clusterposx[i]*clusterspacing)
                FinalStringy.append(np.sin(clustertheta[i])*stringposx[j]*spacing + np.cos(clustertheta[i])*stringposy[j]*spacing + clusterposy[i]*clusterspacing)

    mean_x = sum(FinalStringx)/len(FinalStringx)
    mean_y = sum(FinalStringy)/len(FinalStringy)

    for i in range(len(FinalStringx)):
        FinalStringx[i] -= mean_x
        FinalStringy[i] -= mean_y

    sp = 950.0/19.0
    depthlist = [(-450.0+sp*i)*I3Units.meter for i in range(20)]
    depth = np.array(depthlist)

    for i in range(len(FinalStringx)) :
        for m in range(domsPerString):
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(FinalStringx[i], FinalStringy[i], depth[m])
            for j in range(args.npmts) :
                omkey = OMKey(i+1, m+1, j+1)
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
