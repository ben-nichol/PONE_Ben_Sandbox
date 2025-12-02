from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--spacing", type=float, default=80.0, help="Spacing for strings in meters.")
parser.add_argument("-d", "--ndoms", type=int, default=20, help="DOMs per string.")
parser.add_argument("-r", "--domradius", type=float, default=(17.0 * 2.54 * 0.01 * 0.5), help='Radius of DOM in meters. Defaults to 17"')
parser.add_argument("-p", "--npmts", type=int, default=16, help="PMTs per DOM.")
parser.add_argument("-o", "--output", type=str, default="Special6String", help="Output filename prefix.")

args = parser.parse_args()

outfileName = f"PONE_{args.output}_{int(args.spacing)}Spacing.i3.gz"
outfile = dataio.I3File(outfileName, "w")

# Generate 6-string geometry positions (central string + 5 surrounding strings)
offset = 0.0
anglediff = np.pi * (2.0 / 5)
neighbourangles = [offset + i * anglediff for i in range(5)]

xpositions = [0.0] + [np.cos(angle) for angle in neighbourangles]
ypositions = [0.0] + [np.sin(angle) for angle in neighbourangles]

# Create depth list
sp = 950.0 / 19.0  # spacing between DOMs
depthlist = [(-450.0 + sp * i) * I3Units.meter for i in range(args.ndoms)]

# Scale positions by spacing
scaled_xpositions = [x * args.spacing for x in xpositions]
scaled_ypositions = [y * args.spacing for y in ypositions]

# Generate frames using gcdHelpers
omsequence = ['mDOM'] * args.ndoms  # All mDOMs for 6-string config
gframe = gcdHelpers.generateGFrame(scaled_xpositions, scaled_ypositions, depthlist, omsequence, args.domradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
