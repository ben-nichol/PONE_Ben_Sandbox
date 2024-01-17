# ---------------------------------------------
# Classes and functions for the K40 noise generator
#
# Written by Jakub Stacho
# Last modified: Nov 23, 2022
# ---------------------------------------------

import numpy as np
import pickle as pickle

from icecube import icetray, dataio, dataclasses, simclasses  # , clsim
from icecube.icetray import I3Units


# -----------------------------------------------------------------------------------------------
def unpickle(filePath):
    """
    Unpickle the file at the given file path
    This is used to unpickle the noise characterization class

    Parameters:
        filePath  - string path to the pickled file
    """
    with open(filePath, "rb") as inp:
        # return(pickle.load(inp))
        return pickle.load(inp, encoding="latin1")


# -----------------------------------------------------------------------------------------------
def makePhoton(position, direction, time, wavelength, minorID):
    """
    Returns an I3Photon with the passed attributes

    Parameters:
    position    - (np.array) - xyz coordinates for photon hit position
    direction   - (np.array) - xyz coordinates for photon incident direction
    time        - (float) - Hit timestamp
    wavelength  - (float) - Wavelength of the incident photon
    minorID     - (int) - Minor ID of the particle
    """

    photon = simclasses.I3Photon()
    # photon = simclasses.I3CompressedPhoton()

    # particle IDs included for coincidence analysis

    photon.time = time
    photon.wavelength = wavelength * 1e-9
    photon.weight = 1
    photon.particleMajorID = 0
    photon.particleMinorID = minorID

    pPos = position * I3Units.m
    photon.pos = dataclasses.I3Position(pPos[0], pPos[1], pPos[2])
    pDir = direction * I3Units.m
    photon.dir = dataclasses.I3Direction(pDir[0], pDir[1], pDir[2])

    return photon


# -----------------------------------------------------------------------------------------------
def getMeanEventRate(s=3.482, frameLen=10000):
    """
    Returns the mean event rate of K40 detections

    Parameters:
        frameLen  - Length of frames used in simulation [ns]
                    Only change if there has been a new
                    K40 noise characterization
        s         - Ocean salinity [%]


    Comments:

    The simulation that this noise characterization is based on was only
    run at a salinity of (3.482%). The salinity impacts the number of
    events seen in each frame which is poisson distributed with a mean of
    25.4 events per 10 us frame. (!!!This is using the KM3NeT DOM!!!)

    I think if we want to have salinity as an input parameter to the noise
    generator we might be able to get away with adjusting this mean.

    Based on the current simulation input parameters, the activity of the
    K40 in the Cascadia Basin is 12.1223 decays /ms /m3.
    (!!! The above might be closer to 12.1290 !!!)

    This corresponds to 0.121223 decays /10us /m3

    For the 50 m radius spherical world volume used,
    this is a total of 63472.214 decays /10us

    Our DOM on average detects 25.4 decays /10us which we can say
    corresponds to an efficiency of 0.04 %

    Hence, if we want to adjust the salinity, we can recalculate the
    total number of decays over some time and thus the mean we
    should expect get detected.
    """
    # ---------------------------------------------
    # Define constants
    rk = 1.11  # Potassium fraction in sea salt [%]
    ri = 0.0117  # K40 isotope fraction [%]
    p = 1.013  # Ocean water density [g/cm3]
    NA = 6.022e23  # Avogadro's Number [/mol]
    A = 39.96  # K40 atomic weight [g/mol]
    hl = 1.251e9  # K40 half life [years]

    R = 50  # Simulation world radius [m]
    e = 0.04  # Detection efficiency [%] - Determined from characterization

    # ---------------------------------------------
    # Calculate
    num = (s / 100) * (rk / 100) * (ri / 100) * p * (1e6) * NA * np.log(2)
    den = A * hl * (365 * 24 * 60 * 60 * 1000)

    Bq = (num / den) / 100  # Convert to /10us

    totFrameEvents = Bq * (4.0 / 3.0) * np.pi * (R**3)
    detFrameEvents = totFrameEvents * (e / 100)
    eventRate = detFrameEvents / frameLen

    return eventRate


# -----------------------------------------------------------------------------------------------
class myQuaternion:
    """
    My own simple quaternion class for rotations
    """

    def __init__(self, axis, angle):
        """
        Creates a quaternion object with the specified rotation axis and angle

        Parameters:
            axis   - np.array xyz components of the rotation axis
            angle  - angle to rotate by in radians
        """
        ang = angle / 2

        # Normalize the rotation axis
        axis /= np.linalg.norm(axis)

        uxyz = axis * np.sin(ang)
        s = np.cos(ang)

        # Normalize
        q = np.array([uxyz[0], uxyz[1], uxyz[2], s])
        q /= np.linalg.norm(q)

        self.u = np.array([q[0], q[1], q[2]])
        self.s = q[3]

    def rotate(self, v):
        """
        Applies this quaternion on a vector and returns the rotated vector

        Parameters:
            v - np.array of a vector that should be rotated
        """

        a = 2 * self.u.dot(v) * self.u
        b = (self.s * self.s - self.u.dot(self.u)) * v
        c = 2 * self.s * np.cross(self.u, v)

        vPrime = a + b + c

        return vPrime

    def printQ(self):
        print(
            "%.3f + %.3fi + %.3fj + %.3fk" % (self.s, self.u[0], self.u[1], self.u[2])
        )


# -----------------------------------------------------------------------------------------------
class K40:
    """
    Class for holding and sampling K40 distributions
    Binnings for distributions are hard coded here whereas the distributions
    themselves are passed in when the class is created
    """

    # Define binnings
    dBs = np.linspace(0, 50, 51)  # Distance Bins
    dBCentres = np.linspace(0.5, 49.5, 50)

    linLBins = np.linspace(0, 1, 51)
    linLBinCentres = np.linspace(0.01, 0.99, 50)
    linAngBins = np.linspace(0, 180, 46)
    linAngBinCentres = np.linspace(2, 178, 45)
    linStdvBins = np.linspace(0, 0.25, 26)
    linStdvBinCentres = np.linspace(0.005, 0.245, 25)

    tDBs = np.linspace(0, 60, 61)  # Hit Time difference Bins
    tDBCentres = np.linspace(0.5, 59.5, 60)

    wlBins = np.linspace(200, 700, 51)  # Wavelength Bins
    wlBinCentres = np.linspace(205, 695, 50)

    angDevBins = np.linspace(0, 180, 91)  # Angular Deviation Bins
    angDevBinCentres = np.linspace(0.5, 179.5, 90)

    singleIncAngBins = np.linspace(0, 90, 46)  # Incident Angle Bins for Single Hits
    singleIncAngBinCentres = np.linspace(0.5, 89.5, 45)

    bigHitBins = np.linspace(6.5, 46.5, 41)  # Large n Photon Event Photon Number Bins
    bigHitBinCentres = np.linspace(7, 46, 40)

    # Binnings for distributions that depend on zenith and distance to decay
    zenDZBins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    zenDistBins = [0, 2, 4, 6, 10, 50]
    angDevZBins = [0, 30, 60, 90]

    def __init__(self, dD, sP, nP, bN, sD, pD, eAD, wls, aDv, sID, tDiff):
        """
        Initialize the class with the K40 characterization
        """
        self.dDist = dD  # Total (all # all scatter) distance to decay probability distribution [50]
        self.sProbs = sP  # Scatter probabilities as a function of distance [50][4]
        self.nPProbs = (
            nP  # Number of photon probability as a function of distance [50][10]
        )
        self.bigHitN = (
            bN  # Number of photons distributions for events with >6 photons [4][40]
        )
        self.linSpreadDist = sD  # Point Spread Transverse [6].... ??????
        self.linPointDist = pD  # Point Spread Longitudinal [3][50]
        self.eventAngDist = eAD  # Angle between circle endpoints [4][45]
        self.wlDists = wls  # Wavelength distribution [5][50]
        self.angDevDist = aDv  # Incident angle deviation distribution [3][3][90]
        self.singleIncDist = sID  # Incident angle distribution for single hits [45]
        self.tDiffs = (
            tDiff  # Arrival time difference distribution [4][...] for each scattering
        )
        # self.eRate = eR          # Event rate - events/ns - No longer fixed salinity

    def getBin(self, val, bins):
        """
        Returns the number(index) of the bin that the given value falls into
        """
        binNum = np.digitize(val, bins)
        # Need to also include the last bin edge manually
        if val >= bins[-1]:
            binNum = binNum - 1
        return binNum - 1

    def smear(self, val, bins):
        """
        Returns a given value smeared by a random number within half the distribution bin size
        """
        s = abs(bins[0] - bins[1]) / 2
        shift = np.random.uniform(-s, s)
        return val + shift

    def getDir(self, num, dim=3):
        """
        Generates an array of unit vector in a random direction for a given number of spatial dimensions
        """
        dirs = np.array([None for i in np.arange(num)])
        for i in np.arange(num):
            p = np.random.randn(dim)
            while np.linalg.norm(p) < 0.01:
                p = np.random.randn(dim)
            p /= np.linalg.norm(p)
            dirs[i] = p

        return dirs

    def getD(self, num):
        """
        Returns an array of distances to a K40 decay and the corresponding index
        """
        ds = np.array([None for i in np.arange(num)])
        for i in np.arange(num):
            dIndex = np.random.choice(np.arange(len(self.dBCentres)), 1, p=self.dDist)[
                0
            ]
            dSampled = self.dBCentres[dIndex]
            d = self.smear(dSampled, self.dBCentres)
            ds[i] = (dIndex, d)

        return ds

    def getS(self, dIndex):  # NO LONGER USED TO AVOID FUNCTION CALL OVERHEAD
        """
        Returns the number of scatters based on distance to decay
        """
        ss = np.random.choice(np.arange(4), 1, p=self.sProbs[dIndex])[0]
        return ss

    def getNumPs(self, ds, num):
        """
        Returns an array of the number of photons based on distance to decay
        """
        nPArr = np.arange(num)

        for i in np.arange(num):
            nPs = np.random.choice(np.arange(10), 1, p=self.nPProbs[ds[i][0]])[0] + 1

            # If one of the big hit bins is selected we need to sample the big hit distribution
            # to determine the final number of photons that hit the DOM
            if nPs > 6:
                nPs = int(
                    np.random.choice(self.bigHitBinCentres, 1, p=self.bigHitN[nPs - 7])[
                        0
                    ]
                )
            nPArr[i] = nPs

        return nPArr

    def getLinPoints(self, nPs):
        """
        Returns the start and end points of the hit line on a projected 2d circle
        as a function of the number of photons in the event
        """
        # Choose a random point on the circle
        s = self.getDir(num=1, dim=2)[0]

        # Choose the bin based on # of photons in the event
        if nPs < 10:
            ind = 0
        elif nPs < 15:
            ind = 1
        elif nPs < 25:
            ind = 2
        else:
            ind = 3

        # Get second point
        ang = self.smear(
            np.random.choice(self.linAngBinCentres, 1, p=self.eventAngDist[ind])[0],
            self.linAngBinCentres,
        )
        ang = ang * np.pi / 180
        eX = np.cos(ang) * s[0] - np.sin(ang) * s[1]
        eY = np.cos(ang) * s[1] + np.sin(ang) * s[0]
        e = np.array([eX, eY])

        return s, e

    def getOffsetRIndex(self, start, end):
        """
        Returns an offset index based on the closest radial approach from a line to the
        centre of the circle.
        """
        # Get closest approach to centre
        o = np.array([0, 0])
        originOffset = np.linalg.norm(
            np.cross(end - start, start - o)
        ) / np.linalg.norm(end - start)

        offsetIndex = 2
        if originOffset < 1.0 / 3.0:
            offsetIndex = 0
        elif originOffset < 2.0 / 3.0:
            offsetIndex = 1

        return offsetIndex

    def getOffsetDIndex(self, d2d):
        """
        Returns an offset index based on the distance to decay used for generating
        points along a line.
        """
        if d2d > 5.0:
            offsetIndex = 5
        elif d2d > 2.0:
            offsetIndex = 4
        elif d2d > 1.0:
            offsetIndex = 3
        elif d2d > 0.5:
            offsetIndex = 2
        elif d2d > 0.2:
            offsetIndex = 1
        else:
            offsetIndex = 0

        return offsetIndex

    def getStdv(self, offsetDIndex):  # NO LONGER USED TO AVOID FUNCTION CALL OVERHEAD
        """
        Returns the transverse standard deviation for points along a line
        given an offsetDIndex reffering to the distance to the decay
        """
        tStdv = np.random.choice(
            self.linStdvBinCentres, 1, p=self.linSpreadDist[offsetDIndex]
        )[0]
        return tStdv

    def getLinPos(self, offsetRIndex):  # NO LONGER USED TO AVOID FUNCTION CALL OVERHEAD
        """
        Returns the normalized position of a hit along a given line given
        an offsetRIndex reffering to the radial distance between the
        line and the centre of the projected circle
        """
        lPos = self.smear(
            np.random.choice(self.linLBinCentres, 1, p=self.linPointDist[offsetRIndex])[
                0
            ],
            self.linLBinCentres,
        )
        return lPos

    def getWl(self, d):
        """
        Returns a photon wavelength based on distance to decay
        """
        dBin = self.getBin(d, self.zenDistBins)
        wlSampled = np.random.choice(self.wlBinCentres, 1, p=self.wlDists[dBin])[0]
        wl = self.smear(wlSampled, self.wlBinCentres)
        return wl

    def getDevAng(self, z, s):
        """
        Returns the deviation of incident angle from the direct line for a
        photon with a zenith of <= 90
        """
        if s == 0:
            aDev = 0
        else:
            zBin = self.getBin(z, self.angDevZBins)
            devSampled = np.random.choice(
                self.angDevBinCentres, 1, p=self.angDevDist[zBin][s - 1]
            )[0]
            aDev = self.smear(devSampled, self.angDevBinCentres)
        return aDev

    def getIncDirNear(self, hitPos, dist, angDev):
        """
        Returns the incident direction of a photon with a zenith of < 90 given
        the position of the photon, the distance to the decay, and the incident
        deviation from the direct line
        """
        origin = np.array([0, 0, dist])
        direct = hitPos - origin
        directDir = direct / np.sqrt(
            direct.dot(direct)
        )  # Get the direct line direction

        zenith = np.arccos(directDir[2])
        newZenith = (
            ((zenith * 180 / np.pi) + angDev) * np.pi / 180
        )  # Determine the new zenith

        newZ = np.cos(newZenith)
        xyScale = np.sin(newZenith) / np.sin(
            zenith
        )  # Scale x and y to z so that norm is 1
        newX = directDir[0] * xyScale
        newY = directDir[1] * xyScale

        return np.array([newX, newY, newZ])

    def getSingleInc(self):
        """
        Returns a unit vector representing the incident angle for a
        single photon hit
        """
        # Sample an incident angle in degrees
        ang = self.smear(
            np.random.choice(self.singleIncAngBinCentres, 1, p=self.singleIncDist)[0],
            self.singleIncAngBinCentres,
        )
        ang = ang * np.pi / 180

        # Rotate towards x by the given angle
        inc = np.array([np.sin(ang), 0, -np.cos(ang)])

        # Choose another angle and rotate in azimuth
        az = np.random.uniform(0, 2 * np.pi)
        incDir = np.array(
            [
                inc[0] * np.cos(az) - inc[1] * np.sin(az),
                inc[0] * np.sin(az) + inc[1] * np.cos(az),
                inc[2],
            ]
        )

        # Normalize again just in case
        incDir /= np.linalg.norm(incDir)

        return incDir

    def getNextTime(self, t, s, lowerT, upperT):
        """
        Returns the subsequent time of a hit based on the time of the
        previous hit as well as # scatters (Framelength in ns)
        """
        tDiffSampled = np.random.choice(self.tDBCentres, 1, p=self.tDiffs[s])[0]
        direction = 1 if np.random.random() < 0.5 else -1
        tDiff = direction * self.smear(tDiffSampled, self.tDBCentres)

        newT = t + tDiff

        # Make sure the photons remain in the frame
        if newT < lowerT:
            newT = lowerT
        if newT > upperT:
            newT = upperT

        return newT

    def getFrameEvents(self, frameLength, eRate):
        """
        Returns a list of event times within a frame of a given length (in ns)
        and an event rate in (events/ns)
        """
        frameEs = []

        # Pick some time for the first event in the frame
        currentTime = np.random.gamma(1, 1 / eRate)

        while currentTime < frameLength:
            frameEs.append(currentTime)

            # Get the time to the next event
            eDiff = np.random.gamma(1, 1 / eRate)

            # Update the current time
            currentTime = currentTime + eDiff

        return np.array(frameEs)
