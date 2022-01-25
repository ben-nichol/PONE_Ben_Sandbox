"""!
DOM Utilities is a collection of functions and variables for the DOMs. 
"""
import numpy as np
import os
from icecube.icetray import I3Units, OMKey
from icecube import icetray, dataclasses, dataio, simclasses
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#List for PMT acceptance table
PMTacceptance = list()
#Overall acceptance
Totalacceptance = list()
#List for PMT directions
PMTDirection = list()
#Maximum value for the PMT acceptance
maxAngularAcceptance = 0.0
#Max value for the quantum efficiency
maxQE = 0.0
#List for PMT QE
PMTQE = list()
#Boolian for if acceptance parameters are loaded
AcceptanceLoaded = False
#Boolian for if PMT QE parameters are loaded
QELoaded = False
#random service
random = None

def SetRandomService(_RandomService):
    global random
    random = _RandomService

def AddCLSimPhotonWeighttoFrame(frame,photonweight=1) :
    frame["CLSimPhotonWeight"] = dataclasses.I3Double(photonweight)
    return True

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
def GetPMTAcceptance(infile = os.getenv('PONESRCDIR')+"/data/PMTAcceptance_13PMTConfig.txt") :
    global PMTacceptance
    global PMTDirection
    global maxAngularAcceptance
    global AcceptanceLoaded

    domaccFile = open(infile,"r")
    lines = domaccFile.readlines()
    maxTotaleff = 0.0;

    PMTacceptance = list()
    PMTDirection = list()
    maxAngularAcceptance = 0.0
    Totalacceptance = list()

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

    for i in range(len(PMTacceptance[0])) :
        Totalacceptance.append([])
        for j in range(len(PMTacceptance[0][0])) :
            Totalacceptance[-1].append(0.0)
            for n in range(len(PMTacceptance)) :
                    Totalacceptance[-1][-1] += PMTacceptance[n][i][j]
            if Totalacceptance[-1][-1] > maxAngularAcceptance :
                maxAngularAcceptance = Totalacceptance[-1][-1]

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
    AcceptanceLoaded = True


"""!
GetPMTDirection(pmtid)
Iputs: pmtid = the PMT number within the DOM (OMKey.pmt for IceTray pulseseries keys)
Operation:
		Returns the x,y,z for the PMT direction. 
"""
def GetPMTDirection(pmtid):

    global PMTDirection
    global AcceptanceLoaded
    if not AcceptanceLoaded:
        GetPMTAcceptance()
    x = np.sin(PMTDirection[pmtid][0])*np.cos(PMTDirection[pmtid][1])
    y = np.sin(PMTDirection[pmtid][0])*np.sin(PMTDirection[pmtid][1])
    z = np.cos(PMTDirection[pmtid][0])
    return x, y, z

"""!
GetPMTQETable(infile)
Inputs: infile = file that contains the PMT QE information. 
Operation:
    Reads the file and buils the PMT QE table.
"""
def GetPMTQETable(infile = os.getenv('PONESRCDIR')+"/data/PMTQE.txt") :
    global PMTQE
    global maxQE
    global QELoaded

    PMTQE_value = list()
    PMTQE_wl = list()
    PMTQE = list()
    maxQE = 0.0

    pmtqeFile = open(infile,"r")
    lines = pmtqeFile.readlines()
    for line in lines :
        splitline = line.split(",",100)
        PMTQE_wl.append(float(splitline[0]))
        PMTQE_value.append(float(splitline[1]))
    j = 1
    #Take the data in the file and make an equal spaced table
    #extrapolating between the points. 
    for i in range(1000):
        if i<PMTQE_wl[0] :
            PMTQE.append(0.0)
            continue
        elif i>PMTQE_wl[-1] :
            PMTQE.append(0.0)
            continue
        
        while i > PMTQE_wl[j] :
            j += 1
        QE = PMTQE_value[j-1] + (PMTQE_value[j]-PMTQE_value[j-1])*((float(i)-PMTQE_wl[j-1])/(PMTQE_wl[j]-PMTQE_wl[j-1]))
        PMTQE.append(QE)
        if QE > maxQE :
            maxQE = QE

    QELoaded = True


def MakePMTAcceptancePlots(directory):

    global PMTacceptance
    global maxAngularAcceptance
    global AcceptanceLoaded
    if not AcceptanceLoaded:
        GetPMTAcceptance()

    print(maxAngularAcceptance)

    for n in range(len(PMTacceptance)) :
        pmtprobs = []
        x = []
        y = []

        for i in range(179) :
            for j in range(359) :
                x.append((1./180.)*np.pi*j)
                y.append((1./180.)*np.pi*i)
                pmtprobs.append(PMTacceptance[n][i][j]/maxAngularAcceptance)

        plt.hist2d(x, y, bins=(358, 178), cmap=plt.cm.jet,weights=pmtprobs)
        plt.savefig(directory+"/PMTAcceptance_"+str(n)+".png")

def MakeQEPlot(directory) :
    global PMTQE
    global maxQE
    global QELoaded
    if not QELoaded :
        GetPMTQETable()
    plt.plot(PMTQE)
    plt.savefig(directory+"/QuantumEfficiency.png")
    
"""!
GetPMTQE(wl)
Inputs:
    wl = wavelength (m)
Operation:
    Returns the quantum efficiency 
"""
def GetPMTQE(wl):
    global PMTQE
    global QELoaded

    if not QELoaded :
        GetPMTQETable()

    if int(wl*1.0e9) > len(PMTQE)-1 :
        return 0.0

    return PMTQE[int(wl*1.0e9)]

"""!
GetPMTQEnm(wl)
Inputs:
    wl = wavelength (nm)
Operation:
    Returns the quantum efficiency, this time taking in a wavelength in nm.
"""
def GetPMTQEnm(wl):
    global PMTQE
    global QELoaded

    if not QELoaded :
        GetPMTQETable()

    if int(wl) > len(PMTQE)-1 :
        return 0.0

    return PMTQE[int(wl)]

"""!
GetMaxTotalAcceptance()
Inputs: non
Operation:
    Returns the maximum total acceptance for the DOM, this helps scale CLSim to make it more efficient. 
"""
def GetMaxTotalAcceptance() :
    global maxQE
    global maxAngularAcceptance
    global AcceptanceLoaded
    global QELoaded

    if not QELoaded :
        GetPMTQETable()

    if not AcceptanceLoaded:
        GetPMTAcceptance()

    return maxQE*maxAngularAcceptance

"""!
GetMaxAngularAcceptance()
Inputs: NONE
Operation:
    Returns the maximum value from the angular acceptance table.
"""

def GetMaxAngularAcceptance() :
    global maxAngularAcceptance
    global AcceptanceLoaded

    if not AcceptanceLoaded:
        GetPMTAcceptance()
    return maxAngularAcceptance

def GetMaxPMTQE() :
    global maxQE
    global QELoaded

    if not QELoaded :
        GetPMTQETable()
    return maxQE

"""!
GetPMT(photonDir,wl,random)
Input: 
    photonDir = list of x,y,z or theta,phi
    wl = wavelength (m or nm)
    random = random number between 0 and 1
Operation:
    Randomly assigns photon hits to PMTs based on the photon direction and PMT acceptance.
    This function also applies the QE shape. This function assumes that CLSim has already
    scaled down the photon production so that the maxacceptance after CLSim is scaled to 1.0.

"""
def GetPMT(photonDir,wl,photonweight):
    global PMTacceptance
    global maxQE
    global maxAngularAcceptance 
    global AcceptanceLoaded

    if not AcceptanceLoaded:
        GetPMTAcceptance()

    QEweight = maxQE/photonweight
    Angweight = 1.0
    if QEweight > 1 :
        Angweight = 1.0/QEweight
        QEweight = 1.0


    QEProb = 0.0
    if wl < 0.1 :
        QEProb = GetPMTQE(wl)/QEweight
    else :
        QEProb = GetPMTQEnm(wl)/QEweight


    if random.uniform(0.0,1.0) > QEProb :
        return -1
  
    theta = 0.0
    phi = 0.0
    if len(photonDir) < 3 :
        theta = photonDir[0]
        phi = photonDir[1]
    else :
        theta = np.arccos(photonDir[2])
        phi = np.arctan2(photonDir[1],photonDir[0])

    thetaBin = max(0,min(178,int(180.0*theta/np.pi)))
    phiBin = max(0,min(358,int(180.0*phi/np.pi)))
    pmtprobs = []
    for i in range(len(PMTacceptance)) :
        pmtprobs.append(PMTacceptance[i][thetaBin][phiBin]/Angweight)

    totalprob = sum(pmtprobs)
    rand = random.uniform(0.0,1.0)
    if rand > totalprob :
        return -1
    i=0;

    sumprob = pmtprobs[i]
    while rand > sumprob/totalprob and i<len(pmtprobs)-1 :
       i+=1
       sumprob += pmtprobs[i]

    return i

def GetNPMTs() :
    global PMTacceptance
    global AcceptanceLoaded

    if not AcceptanceLoaded:
        GetPMTAcceptance()
    return len(PMTacceptance)

def NoPMTKey(domkey) :
    return OMKey(domkey.string,domkey.om,0)

def AddPMTKey(domkey,ipmt) :
    return OMKey(domkey.string,domkey.om,ipmt)
