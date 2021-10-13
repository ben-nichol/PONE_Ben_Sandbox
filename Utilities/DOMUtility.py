"""!
DOM Utilities is a collection of functions and variables for the DOMs. 
"""
import numpy as np

#List for PMT acceptance table
PMTacceptance = list()
#List for PMT directions
PMTDirection = list()

"""!
GetPMTAcceptance(infile)
Inputs: infile = text file with photon acceptance information for the PMTs in a DOM.
Operation:
	The Get PMT acceptance method reads in the text file of photon acceptances for 
	plane waves and builds the PMT acceptance table. The method then finds the centroid 
	in the acceptance direction and assigns this as the PMT's viewing direction. Please 
	note that the PMT direction if defined by the photon travel direction and is opposite
	the PMT mounting direction.
"""
def GetPMTAcceptance(infile) :
	global PMTacceptance
	global PMTDirection

   	domaccFile = open(infile,"r")
    lines = domaccFile.readlines()
    maxTotaleff = 0.0;

    PMTacceptance = list()
    PMTDirection = list()

    zenithcount = 0
    for line in lines :
        splitline = line.split(" ",1000)
        if zenithcount % 179 == 0 :
            zenithcount = 0
            PMTacceptance.append([])
        PMTacceptance[-1].append([])
        for value in splitline :
            PMTacceptance[-1][-1].append(float(value))
        zenithcount += 1

        #Compute PMT view direction from the centroid of the acceptance.
        for i in range(len(PMTacceptance)):
            acceptancesum = 0.0
            x,y,z = 0.0,0.0,0.0
            for theta in range(len(PMTacceptance[i])):
                for phi in range(len(PMTacceptance[i][theta])):
                    x += np.sin(float(theta)*np.pi/180.)*np.cos(float(phi)*np.pi/180.)*PMTacceptance[i][theta][phi]
                    y += np.sin(float(theta)*np.pi/180.)*np.sin(float(phi)*np.pi/180.)*PMTacceptance[i][theta][phi]
                    z += np.cos(float(theta)*np.pi/180.)*PMTacceptance[i][theta][phi]
                    acceptancesum += PMTacceptance[i][theta][phi]
            x /= acceptancesum
            y /= acceptancesum
            z /= acceptancesum
            r = np.sqrt(x*x+y*y+z*z)
            x /= r
            y /= r
            z /= r
            theta = np.arccos(z)
            phi = 0.0
            if theta != 0.0 or theta != np.pi :
                phi = np.arctan2(y,x)

            if phi < 0.0 :
                phi += 2.0*np.pi
            PMTDirection.append([theta,phi])


"""!
GetPMTDirection(pmtid)
Iputs: pmtid = the PMT number within the DOM (OMKey.pmt for IceTray pulseseries keys)
Operation:
		Returns the 
"""
def GetPMTDirection(pmtid):

	global PMTDirection
	x = np.sin(PMTDirection[pmtid][0])*np.cos(PMTDirection[pmtid][1])
	y = np.sin(PMTDirection[pmtid][0])*np.sin(PMTDirection[pmtid][1])
	z = np.cos(PMTDirection[pmtid][0])
    return x, y, z

