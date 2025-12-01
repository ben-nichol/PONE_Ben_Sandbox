from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
from LatticeCalculator import generateLatticeSpots
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--spacing", type=float, default=80.0, help="Spacing for strings within cluster (meters).")
parser.add_argument("-l", "--clusterspacing", type=float, default=400.0, help="Spacing between clusters (meters).")
parser.add_argument("-n", "--nstring", type=int, default=10, help="Number of strings per cluster.")
parser.add_argument("-c", "--nclusters", type=int, default=7, help="Number of clusters.")
parser.add_argument("-d", "--ndoms", type=int, default=20, help="DOMs per string.")
parser.add_argument("-r", "--domradius", type=float, default=(17.0 * 2.54 * 0.01 * 0.5), help='Radius of DOM in meters. Defaults to 17"')
parser.add_argument("-o", "--output", type=str, default=None, help="Output filename prefix.")

args = parser.parse_args()

# Generate output filename
if args.output is None:
    outfileName = f"PONE_{args.nstring}String_{args.nclusters}Cluster.i3.gz"
else:
    outfileName = f"PONE_{args.output}.i3.gz"

outfile = dataio.I3File(outfileName, "w")

# Generate lattice positions for strings within clusters and cluster positions
stringposx, stringposy, theta = generateLatticeSpots(args.nstring)
clusterposx, clusterposy, clustertheta = generateLatticeSpots(args.nclusters)

# Calculate final string positions by combining cluster and string layouts
final_xpositions = []
final_ypositions = []

for i in range(len(clusterposx)):
    for j in range(len(stringposx)):
        if i == 0:  # First cluster - no rotation
            final_x = stringposx[j] * args.spacing + clusterposx[i] * args.clusterspacing
            final_y = stringposy[j] * args.spacing + clusterposy[i] * args.clusterspacing
        else:  # Rotate strings according to cluster orientation
            cos_theta = np.cos(clustertheta[i])
            sin_theta = np.sin(clustertheta[i])
            rotated_x = cos_theta * stringposx[j] * args.spacing - sin_theta * stringposy[j] * args.spacing
            rotated_y = sin_theta * stringposx[j] * args.spacing + cos_theta * stringposy[j] * args.spacing
            final_x = rotated_x + clusterposx[i] * args.clusterspacing
            final_y = rotated_y + clusterposy[i] * args.clusterspacing
        
        final_xpositions.append(final_x)
        final_ypositions.append(final_y)

# Center the layout by removing the mean position
mean_x = sum(final_xpositions) / len(final_xpositions)
mean_y = sum(final_ypositions) / len(final_ypositions)
final_xpositions = [x - mean_x for x in final_xpositions]
final_ypositions = [y - mean_y for y in final_ypositions]

# Create depth list
sp = 950.0 / 19.0  # spacing between DOMs
depthlist = [(-450.0 + sp * i) * I3Units.meter for i in range(args.ndoms)]

# Generate frames using gcdHelpers
omsequence = ['mDOM'] * args.ndoms  # All mDOMs for cluster config
gframe = gcdHelpers.generateGFrame(final_xpositions, final_ypositions, depthlist, omsequence, args.domradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
