from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default=None, help="CSV file to read x and y positions of strings")
parser.add_argument("-o", "--output", type=str, default="out", help="I3 file name to write")
parser.add_argument("-d", "--ndoms", type=int, default=20, help="DOMs per string.")
parser.add_argument("-r", "--domradius", type=float, default=0.2159, help='Radius of DOM in meters. Defaults to 0.2159 m')
parser.add_argument("-l", "--mooringlength", type=int, default=1000, help="Length of the mooring line in meters")
parser.add_argument("-s", "--omsequence", type=list, default=['POM'], help="Repeating sequence of POM or PCAL going from the bottom of the line to the top. Default is all POMs")

args = parser.parse_args()

if args.input is None:
    raise Exception("No input CSV provided. Please provide correct file using -i option")

outfileName = f"PONE_{args.output}.i3.gz"
outfile = dataio.I3File(outfileName, "w")

# Create list of depths for modules
sp = args.mooringlength / args.ndoms  # spacing between DOMs
depthlist = [(sp + sp * i) * I3Units.meter for i in range(args.ndoms)]  # from sp to mooringlength. 0 at sea floor

# Read x,y positions from CSV file
xpositions = []
ypositions = []

with open(args.input, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header if there is one
    for row in reader:
        xpositions.append(float(row[0]))
        ypositions.append(float(row[1]))

# Generate frames using gcdHelpers
gframe = gcdHelpers.generateGFrame(xpositions, ypositions, depthlist, args.omsequence, args.domradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
