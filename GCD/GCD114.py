from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers
import csv
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",  type=str, default="test_pos.csv", help="csv file to read")
parser.add_argument("-o","--output", type=str, default="out", help="i3 file name to write")
parser.add_argument("-d","--noms",  type=int, default=20, help="OMs per string.")
parser.add_argument("-r","--domradius", type=int, default=216, help='Radius of dom. Defaults to 0.216 m')
# parser.add_argument("-p", "--npmts", type=int, default=16, help="PMTs per DOM.")
parser.add_argument("-l", "--mooringlength", type=int, default=1000, help="Length of the mooring")
# parser.add_argument("-s", "--omsequence", type=list, default=['pom'], help="Repeating sequence of pom or pcal going from the bottom of the line to the top. Default is all poms")
parser.add_argument("-s", "--omsequence", type=list, default=['pom','pom','pom','pom','pom','pcal','pom','pom','pom','pom','pom','pom','pcal','pom','pom','pom','pom','pom','pom'], help="Repeating sequence of pom or pcal going from the bottom of the line to the top. Default is all poms")


args = parser.parse_args()

outfileName = (
    "PONE_" + str(args.output)+".i3.gz"
)

outfile = dataio.I3File(outfileName, "w")
domsPerString = args.noms
omsequence = args.omsequence

# MOVE THIS PART TO A MODULE TO LOAD IN. MAYBE INTO POMUtility?
from Utilities.POMUtility import POMProperties
pom_properties = POMProperties()
pmt_vectors = np.array([pom_properties.GetPMTDirection(i+1) for i in range(16)])
PMT_orientations = np.array([dataclasses.I3Orientation(dataclasses.I3Direction(pmt_vectors[i][0],pmt_vectors[i][1],pmt_vectors[i][2])) for i in range(len(pmt_vectors))])



def AddPOM(geomap,position,string,om,separation=0.0785,pom_radius=0.216,npmts=16):
    '''
    This function takes a geomap, xyz position, string number, and om number and adds a POM to that position. 
    The POM is make of two different spheres each corresponding to a different hemisphere. 
    At the moment, this model has each hemisphere facing the x direction and does not have the ability to rotate.

    Parameters
    ----------
    geomap :: dataclasses.I3OMGeoMap() object
        This geomap stores the geometry information which will be modified by this function.
    position :: array like
        position = [x, y, z] XYZ position of the added POM.
    string :: int
        String number.
    om :: int
        OM number.
    separation :: float (default = 0.0785)
        Separation of the centre of the two hemispheres in meters.
    pom_radius :: float (default = 0.216)
        Radius of each hemistphere in meters.
    npmts :: int (default = 16)
        Number of PMTs in the POM.
    '''
    x,y,z = position
    for k in range(2):
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
        # omGeometry.area = area
        omGeometry.area = 4*np.pi*pom_radius**2
        omGeometry.position = dataclasses.I3Position(
            x+separation*(-1)**(k+1), y, z
        )
        for j in range(int(npmts/2)):
            omGeometry.orientation = PMT_orientations[int(j+k*npmts/2)]
            omkey = OMKey(string, om, int(j+k*npmts/2)+1)

            geomap[omkey] = omGeometry

def AddPCAL(geomap,position,string,om,separation=0.0785,pcal_radius=0.216,npmts=8):
    '''
    This function takes a geomap, xyz position, string number, and om number and adds a POM to that position. 
    The POM is make of two different spheres each corresponding to a different hemisphere. 
    At the moment, this model has each hemisphere facing the x direction and does not have the ability to rotate.

    Parameters
    ----------
    geomap :: dataclasses.I3OMGeoMap() object
        This geomap stores the geometry information which will be modified by this function
    position :: array like
        position = [x, y, z] XYZ position of the added PCAL.
    string :: int
        String number.
    om :: int
        OM number.
    separation :: float (default = 0.0785)
        Separation of the two hemispheres in meters.
    pcal_radius :: float (default = 0.216)
        Radius of each hemistphere in meters.
    npmts :: int (default = 8)
        Number of PMTs in the PCAL.
    '''
    warnings.warn("AddPCAL NOT READY YET, NEED TO BE ABLE TO ADD FLASHER TO PROPER LOCATION. HOW DO WE WANT TO NUMBER THE PCAL PMTS?")
    x,y,z = position
    for k in range(2):
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
        # omGeometry.orientation = orientation
        # omGeometry.area = area
        omGeometry.area = 4*np.pi*pcal_radius**2
        omGeometry.position = dataclasses.I3Position(
            x+separation*(-1)**(k+1), y, z
        )
        for j in range(int(npmts/2)):
            omGeometry.orientation = PMT_orientations[int(j+k*8+4)]
            omkey = OMKey(string, om, int(j+k*npmts/2)+5)
            geomap[omkey] = omGeometry

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

    area = 4.0 * ((args.domradius) ** 2.0) * np.pi * I3Units.meter2 
    geomap = dataclasses.I3OMGeoMap()

    #Add OMs at x,y and depth
    for i in range(len(xpositions)):
        for m in range(domsPerString):
            loc = [xpositions[i], ypositions[i], depth[m]]
            string_num = i+1; om_num = m+1
            if omsequence[m%len(omsequence)]=='pom':
                AddPOM(geomap,loc,string_num,om_num)
            elif omsequence[m%len(omsequence)]=='pcal':
                AddPCAL(geomap,loc,string_num,om_num)
            else:
                raise Exception("Unknown optical module type specified. Type can be 'pom' or 'pcal'") 
    
    return geomap


#create list of depths for modules
sp = args.mooringlength / args.noms #spacing
depthlist = [(sp + sp * i) * I3Units.meter for i in range(args.noms)] # from sp to mooringlength. 0 at sea floor
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



orientation = dataclasses.I3Orientation(0, 0, 1, 1, 0, 0) # Making the OM oriented upwards. Don't know if this will change anything

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
ind = 0
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    module = dataclasses.I3ModuleGeo()
    module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
    module.orientation = orientation
    module.pos = dataclasses.I3Position(
            xpositions[dom.string-1], xpositions[dom.string-1], depthlist[dom.om-1]
        )
    # module.radius = np.sqrt(geomap[dom].area / (4.0 * np.pi))
    module.radius = args.domradius
    modgeomap[mkey] = module
    ind+=1
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
