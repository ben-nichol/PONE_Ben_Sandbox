from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import numpy as np
import argparse
import gcdHelpers
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",  type=str, default=None, help="csv file to read in x and y positions of lines")
parser.add_argument("-o","--output", type=str, default="out", help="i3 file name to write")
parser.add_argument("-d","--noms",  type=int, default=20, help="OMs per string.")
parser.add_argument("-r","--pomradius", type=int, default=0.2159, help='Radius of pom. Defaults to 0.2159 m')
parser.add_argument("-l", "--mooringlength", type=int, default=1000, help="Length of the mooring")
parser.add_argument("-s", "--omsequence", type=list, default=['POM'], help="Repeating sequence of POM or PCAL going from the bottom of the line to the top. Default is all POMs")
# The below sequence is the order that PONE-1 will have
# parser.add_argument("-s", "--omsequence", type=list, default=['POM','POM','POM','POM','POM','PCAL','POM','POM','POM','POM','POM','POM','PCAL','POM','POM','POM','POM','POM','POM'], help="Repeating sequence of POM or PCAL going from the bottom of the line to the top. Default is all POMs")


args = parser.parse_args()
if args.input==None:
    raise("No input csv provided. Please provide correct file")

outfileName = (
    "PONE_" + str(args.output)+".i3.gz"
)

outfile = dataio.I3File(outfileName, "w")
domsPerString = args.noms
omsequence = args.omsequence

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


# Generate frames using gcdHelpers
gframe = gcdHelpers.generateGFrame(xpositions, ypositions, depthlist, omsequence, args.pomradius)
geometry = gframe["I3Geometry"]
cframe = gcdHelpers.generateCFrame(geometry, empty=True)
dframe = gcdHelpers.generateDFrame(geometry, empty=True)

outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

outfile.close()
