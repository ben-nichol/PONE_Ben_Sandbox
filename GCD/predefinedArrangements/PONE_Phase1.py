from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default="PONE_Phase1", help="Output filename (without extension)")
parser.add_argument("-c", "--ncircles", type=int, default=2, help="Number of circles in the geometry")
parser.add_argument("-d", "--domsperstring", type=int, default=20, help="Number of DOMs per string")
parser.add_argument("-s", "--stringspercircle", type=int, nargs='+', default=[7, 3], help="Number of strings per circle")
parser.add_argument("-p", "--npmts", type=int, default=13, help="Number of PMTs per mDOM")
parser.add_argument("-r", "--radius", type=float, default=200, help="Maximum radius for circles in meters")

args = parser.parse_args()

outfileName = args.output + ".i3.gz"
outfile = dataio.I3File(outfileName, "w")


def generateGeometry(nCircles, DPS, strings, npmts, max_radius):
    """Generate geometry for PONE Phase 1 detector configuration.
    
    Args:
        nCircles (int): Number of concentric circles
        DPS (int): DOMs per string
        strings (list): Number of strings per circle
        npmts (int): Number of PMTs per mDOM
        max_radius (float): Maximum radius for outermost circle in meters
    
    Returns:
        I3OMGeoMap: Geometry map for all optical modules
    """
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)
    # Area calculation for 17-inch mDOM (diameter = 17 inches = 43.18 cm)
    area = 4.0 * ((17.0 * 2.54 * 0.01 * 0.5) ** 2.0) * np.pi * I3Units.meter2
    geomap = dataclasses.I3OMGeoMap()

    # Generate radii for concentric circles
    radius = np.arange(max_radius, 0, -(max_radius / nCircles))
    
    # Create depth list for DOMs - evenly spaced over 800m mooring
    sp = 800.0 / (DPS - 1)  # Spacing between DOMs
    depthlist = [(-400.0 + sp * i) * I3Units.meter for i in range(DPS)]
    depth = np.array(depthlist)
    
    # Generate x,y positions for all strings
    xPos = []
    yPos = []
    for i in range(len(radius)):
        # Calculate angular spacing for strings on this circle
        spacing = (2 * np.pi * radius[i]) / strings[i]
        thetaDiff = spacing / radius[i]
        theta = np.arange(0, 2 * np.pi, thetaDiff)
        x = radius[i] * np.cos(theta) * I3Units.meter
        y = radius[i] * np.sin(theta) * I3Units.meter
        xPos = np.append(xPos, x)
        yPos = np.append(yPos, y)

    # Create OM geometry for each DOM at each depth
    array = np.ones(xPos.shape)
    zPos = [array * depth[j] for j in range(DPS)]

    for m in range(DPS):  # For each depth level
        for n in range(len(xPos)):  # For each string position
            omGeometry = dataclasses.I3OMGeo()
            omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
            omGeometry.orientation = orientation
            omGeometry.area = area
            omGeometry.position = dataclasses.I3Position(xPos[n], yPos[n], zPos[m][n])
            
            # Create entries for each PMT in the mDOM
            for i in range(npmts):
                omkey = OMKey(n + 1, m + 1, i + 1)
                geomap[omkey] = omGeometry

    return geomap


# Generate geometry
geomap = generateGeometry(args.ncircles, args.domsperstring, args.stringspercircle, args.npmts, args.radius)

# Create geometry object
geometry = dataclasses.I3Geometry()
geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = geomap

# Create G-frame and populate with geometry information
gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
gframe["I3Geometry"] = geometry
gframe["I3OMGeoMap"] = geomap

# Create module geometry map for mDOMs
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

# Create subdetector map
subdetec = dataclasses.I3MapModuleKeyString()
for dom in geomap.keys():
    mkey = dataclasses.ModuleKey(dom.string, dom.om)
    subdetec[mkey] = "PONE"

gframe["Subdetectors"] = subdetec
gframe["StartTime"] = gcdHelpers.start_time
gframe["EndTime"] = gcdHelpers.end_time

# Generate C-frame and D-frame using gcdHelpers
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
