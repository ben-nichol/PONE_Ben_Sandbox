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
parser.add_argument("-l", "--mooringlength", type=int, default=1000, help="Length of the mooring")
parser.add_argument("-s", "--omsequence", type=list, default=['POM'], help="Repeating sequence of POM or PCAL going from the bottom of the line to the top. Default is all POMs")
args = parser.parse_args()


outfileName = (
    "PONE_" + str(args.output)+".i3.gz"
)

outfile = dataio.I3File(outfileName, "w")
domsPerString = args.ndoms
omsequence = args.omsequence

#create list of depths for modules
sp = args.mooringlength / args.ndoms #spacing
depthlist = [(sp + sp * i) * I3Units.meter for i in range(args.ndoms)] # from sp to mooringlength. 0 at sea floor
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

    #Add OMs at x,y and depth
    for i in range(len(xpositions)):
        for m in range(domsPerString):
            loc = [xpositions[i], ypositions[i], depth[m]]
            string_num = i+1; om_num = m+1
            if omsequence[m%len(omsequence)]=='POM':
                gcdHelpers.AddPOM(geomap,loc,string_num,om_num)
            elif omsequence[m%len(omsequence)]=='PCAL':
                gcdHelpers.AddPCAL(geomap,loc,string_num,om_num)
            else:
                raise Exception("Unknown optical module type specified. Type can be 'POM' or 'PCAL'") 
    
    return geomap



geometry = dataclasses.I3Geometry()

geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geomap = generateGeometry()
geometry.omgeo = geomap

gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
gframe["I3Geometry"] = geometry
gframe["I3OMGeoMap"] = geomap

modgeomap = dataclasses.I3ModuleGeoMap()
orientation = dataclasses.I3Orientation(0, 0, 1, 1, 0, 0) # Making the OM oriented upwards. Don't know if this will change anything
ind = 0
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    module = dataclasses.I3ModuleGeo()
    module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
    module.orientation = orientation
    module.pos = dataclasses.I3Position(
            xpositions[dom.string-1], ypositions[dom.string-1], depthlist[dom.om-1]
        )
    module.radius = args.domradius
    modgeomap[mkey] = module
    ind+=1

gframe["I3ModuleGeoMap"] = modgeomap
# Setting subdetector labels to be the POM and PCAL sequence specified
subdetec = dataclasses.I3MapModuleKeyString()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    subdetec[mkey] =omsequence[dom.om%len(omsequence)]
gframe["Subdetectors"] = subdetec
gframe["StartTime"] = gcdHelpers.start_time
gframe["EndTime"] = gcdHelpers.end_time

cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
