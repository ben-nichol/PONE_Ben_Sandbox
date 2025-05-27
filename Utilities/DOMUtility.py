"""!
DOM Utilities is a collection of functions and variables for the DOMs. 
"""
import numpy as np
import os
from icecube.icetray import I3Units, OMKey
from icecube import icetray, dataclasses, dataio, simclasses
from icecube import clsim
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

def NoPMTKey(domkey):
    return OMKey(domkey.string,domkey.om,0)

def AddPMTKey(domkey,ipmt):
    return OMKey(domkey.string,domkey.om,ipmt)

class DOMProperties():
    def __init__(self, 
        PMTAcceptanceFile = os.getenv('PONESRCDIR')+"/data/ang_acceptance_16_all.npy",
        PMTQEFile = os.getenv('PONESRCDIR')+"/data/PMTQE.txt"
        ):

        #List for PMT acceptance table
        self.PMTacceptance = list()

        self.PMTMaxacceptance = list()

        #Overall acceptance
        self.Totalacceptance = list()

        #List for PMT directions
        self.PMTDirection = list()

        #Maximum value for the PMT acceptance
        self.maxAngularAcceptance = 0.0

        #average value of the PMT acceptance
        self.averageAngularAcceptance = 0.0

        #Max value for the quantum efficiency
        self.maxQE = 0.0

        #List for PMT QE
        self.PMTQE = list()

        self._LoadPMTAcceptance(PMTAcceptanceFile)
        self._LoadPMTQETable(PMTQEFile)

    """!
    _LoadPMTAcceptance(infile)
    Inputs: infile = text file with photon acceptance information for the PMTs in a DOM.
    Operation:
        The Get PMT acceptance method reads in the text file of photon acceptances for 
        plane waves and builds the PMT acceptance table. The method then finds the centroid 
        in the acceptance direction and assigns this as the PMT's viewing direction. Please 
        note that the PMT direction if defined by the photon travel direction and is opposite
        the PMT mounting direction.
    """
    def _LoadPMTAcceptance_txt(self, infile) :
        domaccFile = open(infile,"r")
        lines = domaccFile.readlines()
        maxTotaleff = 0.0;

        self.PMTacceptance = list()
        self.PMTMaxacceptance = list()
        self.PMTDirection = list()
        self.maxAngularAcceptance = 0.0
        self.averageAngularAcceptance = 0.0
        self.Totalacceptance = list()

        zenithcount = 0
        for line in lines :
            splitline = line.split(" ",1000)
            if zenithcount % 179 == 0 :
                zenithcount = 0
                self.PMTacceptance.append([])
                self.PMTMaxacceptance.append(0.0)
            self.PMTacceptance[-1].append([])
            for value in splitline :
                self.PMTacceptance[-1][-1].append(float(value))
                self.PMTMaxacceptance[-1] = max(float(value),self.PMTMaxacceptance[-1])
                
            zenithcount += 1

        for i in range(len(self.PMTacceptance[0])) :
            self.Totalacceptance.append([])
            for j in range(len(self.PMTacceptance[0][0])) :
                self.Totalacceptance[-1].append(0.0)
                for n in range(len(self.PMTacceptance)) :
                    self.Totalacceptance[-1][-1] += self.PMTacceptance[n][i][j]
                if self.Totalacceptance[-1][-1] > self.maxAngularAcceptance :
                    self.maxAngularAcceptance = self.Totalacceptance[-1][-1]
        sumtotal = 0.0
        nsum = 0
        #print("max = " + str(max(self.Totalacceptance)))
        self.maxAngularAcceptance = max(self.Totalacceptance)
        for i in range(len(self.Totalacceptance)) :
            sumtotal += sum(self.Totalacceptance[i])
            nsum += len(self.Totalacceptance[i])

        self.averageAngularAcceptance = sumtotal/nsum
        #Compute PMT view direction from the centroid of the acceptance.
        for i in range(len(self.PMTacceptance)):
            acceptancesum = 0.0
            x,y,z = 0.0,0.0,0.0
            for theta in range(len(self.PMTacceptance[i])):
                for phi in range(len(self.PMTacceptance[i][theta])):
                    x += np.sin(float(theta)*np.pi/180.)*np.cos(float(phi)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    y += np.sin(float(theta)*np.pi/180.)*np.sin(float(phi)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    z += np.cos(float(theta)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    acceptancesum += self.PMTacceptance[i][theta][phi]
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
            self.PMTDirection.append([theta,phi])

    def _LoadPMTAcceptance_npy(self, infile) :
        self.PMTacceptance = np.load(infile)

        maxTotaleff = 0.0;

        self.PMTDirection = list()
        self.maxAngularAcceptance = 0.0
        self.Totalacceptance = list()

        for i in range(len(self.PMTacceptance[0])) :
            self.Totalacceptance.append([])
            for j in range(len(self.PMTacceptance[0][0])) :
                self.Totalacceptance[-1].append(0.0)
                for n in range(len(self.PMTacceptance)) :
                    self.Totalacceptance[-1][-1] += self.PMTacceptance[n][i][j]
                if self.Totalacceptance[-1][-1] > self.maxAngularAcceptance :
                    self.maxAngularAcceptance = self.Totalacceptance[-1][-1]

        sumtotal = 0.0
        nsum = 0
        for i in range(len(self.Totalacceptance)) :
            sumtotal += sum(self.Totalacceptance[i])
            nsum += len(self.Totalacceptance[i])
        self.averageAngularAcceptance = sumtotal/nsum

        #Compute PMT view direction from the centroid of the acceptance.
        for i in range(len(self.PMTacceptance)):
            acceptancesum = 0.0
            x,y,z = 0.0,0.0,0.0
            for theta in range(len(self.PMTacceptance[i])):
                for phi in range(len(self.PMTacceptance[i][theta])):
                    x += np.sin(float(theta)*np.pi/180.)*np.cos(float(phi)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    y += np.sin(float(theta)*np.pi/180.)*np.sin(float(phi)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    z += np.cos(float(theta)*np.pi/180.)*self.PMTacceptance[i][theta][phi]
                    acceptancesum += self.PMTacceptance[i][theta][phi]
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
            self.PMTDirection.append([theta,phi])
    

    def _LoadPMTAcceptance(self,infile) :
        extension = infile.split(".",100)[-1]
        if extension == "npy" :
            self._LoadPMTAcceptance_npy(infile)
        else :
            self._LoadPMTAcceptance_txt(infile)



    """!
    GetPMTDirection(pmtid)
    Iputs: pmtid = the PMT number within the DOM (OMKey.pmt for IceTray pulseseries keys)
    Operation:
            Returns the x,y,z for the PMT direction. 
    """
    def GetPMTDirection(self, pmtid):
        x = np.sin(self.PMTDirection[pmtid-1][0])*np.cos(self.PMTDirection[pmtid-1][1])
        y = np.sin(self.PMTDirection[pmtid-1][0])*np.sin(self.PMTDirection[pmtid-1][1])
        z = np.cos(self.PMTDirection[pmtid-1][0])
        return x, y, z

    """!
    _LoadPMTQETable(infile)
    Inputs: infile = file that contains the PMT QE information. 
    Operation:
        Reads the file and buils the PMT QE table.
    """
    def _LoadPMTQETable(self, infile = os.getenv('PONESRCDIR')+"/data/PMTQE.txt") :
        PMTQE_value = list()
        PMTQE_wl = list()
        self.PMTQE = list()
        self.maxQE = 0.0

        pmtqeFile = open(infile,"r")
        lines = pmtqeFile.readlines()
        for line in lines :
            splitline = line.split(",",100)
            PMTQE_wl.append(float(splitline[0]))
            PMTQE_value.append(float(splitline[1]))

        self.wlen_min = np.min(PMTQE_wl)*I3Units.nanometer
        self.wlen_max = np.max(PMTQE_wl)*I3Units.nanometer

        j = 1
        #Take the data in the file and make an equal spaced table
        #extrapolating between the points. 
        for i in range(1000):
            if i<PMTQE_wl[0] :
                self.PMTQE.append(0.0)
                continue
            elif i>PMTQE_wl[-1] :
                self.PMTQE.append(0.0)
                continue
            
            while i > PMTQE_wl[j] :
                j += 1
            QE = PMTQE_value[j-1] + (PMTQE_value[j]-PMTQE_value[j-1])*((float(i)-PMTQE_wl[j-1])/(PMTQE_wl[j]-PMTQE_wl[j-1]))
            self.PMTQE.append(QE)
            if QE > self.maxQE :
                self.maxQE = QE

        self.nominal_clsim_table = self.GetCLSimQETable(factor=1.0)

    def GetCLSimQETable(self, binning = 3.*I3Units.nanometer, factor=1.0):
        num_bins = int((self.wlen_max - self.wlen_min)/binning)

        table = []
        for i in range(num_bins+1):
            wlen = float(i)*binning + self.wlen_min

            lookup_bin = int(wlen/I3Units.nanometer)
            if lookup_bin > len(self.PMTQE)-1:
                value = 0.
            else:
                value = self.PMTQE[lookup_bin]

            table.append(value * factor)

        clsim_table = simclasses.I3CLSimFunctionFromTable(self.wlen_min, binning, table)

        return clsim_table

    def MakePMTAcceptancePlots(self, directory):
        print(self.maxAngularAcceptance)

        for n in range(len(self.PMTacceptance)) :
            pmtprobs = []
            x = []
            y = []

            for i in range(179) :
                for j in range(359) :
                    x.append((1./180.)*np.pi*j)
                    y.append((1./180.)*np.pi*i)
                    pmtprobs.append(self.PMTacceptance[n][i][j]/self.maxAngularAcceptance)

            plt.hist2d(x, y, bins=(358, 178), cmap=plt.cm.jet,weights=pmtprobs)
            plt.savefig(directory+"/PMTAcceptance_"+str(n)+".png")

    def MakeQEPlot(self, directory) :
        plt.plot(self.PMTQE)
        plt.savefig(directory+"/QuantumEfficiency.png")
    
    """!
    GetPMTQE(wl)
    Inputs:
        wl = wavelength (m)
    Operation:
        Returns the quantum efficiency 
    """
    def GetPMTQE(self, wl):
        return self.nominal_clsim_table.GetValue(wl*I3Units.meter)
        # if int(wl*1.0e9) > len(self.PMTQE)-1 :
        #     return 0.0

        # return self.PMTQE[int(wl*1.0e9)]

    """!
    GetPMTQEnm(wl)
    Inputs:
        wl = wavelength (nm)
    Operation:
        Returns the quantum efficiency, this time taking in a wavelength in nm.
    """
    def GetPMTQEnm(self, wl):
        return self.nominal_clsim_table.GetValue(wl*I3Units.nanometer)
        # if int(wl) > len(self.PMTQE)-1 :
        #     return 0.0

        # return self.PMTQE[int(wl)]

    """!
    GetMaxTotalAcceptance()
    Inputs: non
    Operation:
        Returns the maximum total acceptance for the DOM, this helps scale CLSim to make it more efficient. 
    """
    def GetMaxTotalAcceptance(self) :
        return self.maxQE*self.maxAngularAcceptance

    def GetAverageTotalAcceptance(self) :
        return self.averageAngularAcceptance

    """!
    GetMaxAngularAcceptance()
    Inputs: NONE
    Operation:
        Returns the maximum value from the angular acceptance table.
    """
    def GetMaxAngularAcceptance(self) :
        return self.maxAngularAcceptance

    def GetMaxPMTQE(self) :
        return self.maxQE

    def GetNPMTs(self):
        return len(self.PMTacceptance)

    def GetPMTScaledAcceptance(self,pmt,theta,phi) :
        while theta < 0.0 :
            theta += 2.*np.pi
        while theta > 2.0*np.pi :
            theta -= 2.0*np.pi
        if theta > np.pi :
            theta = np.pi - (theta - np.pi)
            phi += np.pi
        while phi < 0.0 :
            phi == 2.*np.pi
        while phi > 2.0*np.pi :
            phi -= 2.0*np.pi

        i = min(max(0,int(theta*180./np.pi)),len(self.PMTacceptance[int(pmt)-1])-1)
        j = min(max(0,int(phi*180./np.pi)),len(self.PMTacceptance[int(pmt)-1][0])-1)
        #print(pmt)
        #print("theta = "+str(theta)+" "+str(i)+" "+str(len(self.PMTacceptance[int(pmt)-1])))
        #print("phi = "+str(phi)+" "+str(j)+" "+str(len(self.PMTacceptance[int(pmt)-1][i])))

        return self.PMTacceptance[int(pmt)-1][i][j]/self.maxAngularAcceptance
