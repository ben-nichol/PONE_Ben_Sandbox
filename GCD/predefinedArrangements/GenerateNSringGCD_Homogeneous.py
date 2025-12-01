from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
from LatticeCalculator import generateLatticeSpots
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--spacing", type=float, default=50.0, help="Spacing for strings (meters).")
parser.add_argument("-n", "--nstring", type=int, default=10, help="Number of strings.")
parser.add_argument("-d", "--ndoms", type=int, default=20, help="DOMs per string.")
parser.add_argument("-r", "--domradius", type=float, default=(17.0 * 2.54 * 0.01 * 0.5), help='Radius of DOM in meters. Defaults to 17"')
parser.add_argument("-o", "--output", type=str, default=None, help="Output filename prefix.")

args = parser.parse_args()

# Generate output filename
if args.output is None:
    outfileName = f"PONE_{args.nstring}String_{int(args.spacing)}Spacing.i3.gz"
else:
    outfileName = f"PONE_{args.output}.i3.gz"

outfile = dataio.I3File(outfileName, "w")

# Generate lattice positions for homogeneous string layout
stringposx, stringposy, theta = generateLatticeSpots(args.nstring)

# Scale positions by spacing
scaled_xpositions = [x * args.spacing for x in stringposx]
scaled_ypositions = [y * args.spacing for y in stringposy]

# Create depth list
sp = 950.0 / 19.0  # spacing between DOMs
depthlist = [(-450.0 + sp * i) * I3Units.meter for i in range(args.ndoms)]

# Optional: Print layout information
print(f"Generated {len(scaled_xpositions)} string positions:")
for i, (x, y) in enumerate(zip(scaled_xpositions, scaled_ypositions)):
    print(f"String {i+1}: ({x:.1f}, {y:.1f}) m")

# Generate frames using gcdHelpers
omsequence = ['mDOM'] * args.ndoms  # All mDOMs for homogeneous config
gframe = gcdHelpers.generateGFrame(scaled_xpositions, scaled_ypositions, depthlist, omsequence, args.domradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
