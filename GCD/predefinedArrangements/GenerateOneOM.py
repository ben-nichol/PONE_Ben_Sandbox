from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--domradius", type=float, default=(17.0 * 2.54 * 0.01 * 0.5), help='Radius of DOM in meters. Defaults to 17"')
parser.add_argument("-o", "--output", type=str, default="one-om-gcd-origin", help="Output filename.")

args = parser.parse_args()

outfileName = f"{args.output}.i3.gz"
outfile = dataio.I3File(outfileName, "w")

# Single string at origin with 3 DOMs
xpositions = [0.0]
ypositions = [0.0]
depthlist = [0.0 * I3Units.meter, 250.0 * I3Units.meter, -250.0 * I3Units.meter]

# Generate frames using gcdHelpers
omsequence = ['mDOM'] * 3  # 3 mDOMs
gframe = gcdHelpers.generateGFrame(xpositions, ypositions, depthlist, omsequence, args.domradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
