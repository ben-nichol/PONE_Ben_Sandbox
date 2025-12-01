#!/usr/bin/env python

"""
a few helper functions that will be commonly used when making
GCD files to test prototype geometries
"""

from AcousticFrame import geometry
from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units
from enum import Enum
import numpy as np
from Utilities.POMUtility import POMProperties
pom_properties = POMProperties()
pmt_vectors = np.array([pom_properties.GetPMTDirection(i+1) for i in range(16)])
PMT_orientations = np.array([dataclasses.I3Orientation(dataclasses.I3Direction(pmt_vectors[i][0],pmt_vectors[i][1],pmt_vectors[i][2])) for i in range(len(pmt_vectors))])

cdfile = dataio.I3File("Calib_and_DetStat_File.i3.gz")
cdframe = cdfile.pop_frame()
calib = cdframe["I3Calibration"]
start_time = calib.start_time
end_time = calib.end_time


def makeCalibrationObject(geometry, empty=False):
    """Creates a calibration with all the correct OMKeys.
    
    Args:
        geometry: The I3Geometry object inputted to the frame.
        empty (bool, optional): If True, returns an empty I3Calibration. Defaults to False.
        
    Returns:
        I3Calibration: An I3Calibration object to be used in the GCD file.
    """
    # # THIS PASSES AN EMPTY I3Calibration FRAME, WORKS FOR ICETRAY V1.14
    if(empty):
        calib = dataclasses.I3Calibration()
        return calib

    # THIS PASSES THE CALIBRATION FRAME TO THE OMKEY(0,0,0) FOR ICETRAY V1.14
    # NEEDS SOME CLEANUP. Why OMKEY(0,0,0)?
    domcal = dataclasses.Map_OMKey_I3DOMCalibration()
    calibrationData = cdframe["I3Calibration"].dom_cal.popitem()[1]
    domcal[OMKey(0,0,0)] = calibrationData
    test = dataclasses.I3Calibration()
    test = cdframe["I3Calibration"]
    test.dom_cal = domcal
    return test


def makeDSObject(geometry, empty=False):
    """Creates a detector status with all the correct OMKeys.
    
    Args:
        geometry: The I3Geometry object inputted to the frame.
        empty (bool, optional): If True, returns an empty I3DetectorStatus. Defaults to False.
        
    Returns:
        I3DetectorStatus: An I3DetectorStatus object to be used in the GCD file.
    """
    if(empty):
        detectorStatus = dataclasses.I3DetectorStatus()
        return detectorStatus
    
    dsData = list(cdframe["I3DetectorStatus"].dom_status.values())[0] # NEW V1.14, NEEDED TO ADD list()
    domstat = dataclasses.Map_OMKey_I3DOMStatus()
    detectorStatus = cdframe["I3DetectorStatus"]
    for omkey in geometry.omgeo.keys():
        domstat[omkey] = dsData

    detectorStatus.dom_status = domstat

    return detectorStatus


def generateCFrame(geometry, empty=False):
    """Creates a frame with all the needed fields in a C frame.
    
    Args:
        geometry: The I3Geometry object inputted to the frame.
        empty (bool, optional): If True, creates frame with empty calibration. Defaults to False.
        
    Returns:
        I3Frame: An I3Frame object to be used as the C frame.
    """
    # intialize frame as a C frame
    frame = icetray.I3Frame(icetray.I3Frame.Calibration)

    # add key-value pairs
    frame["I3Calibration"] = makeCalibrationObject(geometry, empty)
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]

    return frame


def generateDFrame(geometry, empty=False):
    """Creates a frame with all the needed fields in a D frame.
    
    Args:
        geometry: The I3Geometry object inputted to the frame.
        empty (bool, optional): If True, creates frame with empty detector status. Defaults to False.
        
    Returns:
        I3Frame: An I3Frame object to be used as the D frame.
    """
    # intialize frame as a D frame
    frame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

    # add key-value pairs
    frame["I3DetectorStatus"] = makeDSObject(geometry, empty)
    frame["BadDomsList"] = cdframe["BadDomsList"]
    frame["BadDomsListSLC"] = cdframe["BadDomsListSLC"]

    return frame

def makePartialGeometry(denseGeometry, domsUsed):
    """Decomposes a geometry into a new geometry with only the DOMs specified.
    
    Args:
        denseGeometry (I3Geometry): An I3Geometry object representing the dense geometry.
        domsUsed (list): A list of OMKey objects representing the DOMs to be
            included in the new geometry.
            
    Returns:
        I3Geometry: A new I3Geometry object containing only the DOMs in domsUsed.
        
    Raises:
        ValueError: If an OMKey in domsUsed is not in denseGeometry.
    """
    newGeo = dataclasses.I3Geometry()
    newGeo.start_time = start_time
    newGeo.end_time = end_time

    newGeoMap = dataclasses.I3OMGeoMap()
    denseGeoMap = denseGeometry.omgeo
    for omkey in domsUsed:
        if omkey not in denseGeoMap:
            raise ValueError(str(omkey) + "is not in the dense geometry")
        newGeoMap[omkey] = denseGeoMap[omkey]

    newGeo.omgeo = newGeoMap

    return newGeo

def AddPOM(geomap,position,string,om,separation=0.0785,pom_radius=0.2159):
    """Adds a POM to a I3OMGeoMap object at a specified location and IDs.
    
    POM consists of two spheres with correct PMT orientations.
    PMTs 1-8 are on face oriented in the +x direction.
    
    Args:
        geomap (dataclasses.I3OMGeoMap): Object which will have the POM added to it.
        position (array-like): [x,y,z] position of the centre of the POM.
        string (int): String ID that the POM will be added to.
        om (int): OM ID that the POM will be added to.
        separation (float, optional): Separation of the centre of the two 
            hemispheres in meters. Defaults to 0.0785.
        pom_radius (float, optional): Radius of each hemisphere in meters.
            Defaults to 0.2159.
    """
    npmts=16
    x,y,z = position
    for k in range(2):
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
        omGeometry.area = 4*np.pi*pom_radius**2
        omGeometry.position = dataclasses.I3Position(
            x+separation*(-1)**(k), y, z
        )
        for j in range(int(npmts/2)):
            omGeometry.orientation = PMT_orientations[int(j+k*npmts/2)]
            omkey = OMKey(string, om, int(j+k*npmts/2)+1)

            geomap[omkey] = omGeometry

def AddPCAL(geomap,position,string,om,separation=0.0785,pcal_radius=0.2159):
    """Adds a PCAL to a I3OMGeoMap object at a specified location and IDs.
    
    PCAL consists of two spheres with correct PMT orientations.
    PMTs 1-4 are on face oriented in the +x direction.
    
    Args:
        geomap (dataclasses.I3OMGeoMap): Object which will have the PCAL added to it.
        position (array-like): [x,y,z] position of the centre of the PCAL.
        string (int): String ID that the PCAL will be added to.
        om (int): OM ID that the PCAL will be added to.
        separation (float, optional): Separation of the centre of the two 
            hemispheres in meters. Defaults to 0.0785.
        pcal_radius (float, optional): Radius of each hemisphere in meters.
            Defaults to 0.2159.
    """
    npmts=8
    x,y,z = position
    for k in range(2):
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.mDOM
        omGeometry.area = 4*np.pi*pcal_radius**2
        omGeometry.position = dataclasses.I3Position(
            x+separation*(-1)**(k+1), y, z
        )
        for j in range(int(npmts/2)):
            omGeometry.orientation = PMT_orientations[int(j+k*8+4)]
            omkey = OMKey(string, om, int(j+k*npmts/2)+1)
            geomap[omkey] = omGeometry


def generateGeometry(xpositions, ypositions, depth, omsequence):
    """Creates an I3OMGeoMap object from x,y coordinates and other parameters.
    
    Args:
        xpositions (list): List of x positions for strings.
        ypositions (list): List of y positions for strings.
        depth (array): Array of depths for OMs.
        omsequence (list): List defining the OM sequence (POM/PCAL).
        
    Returns:
        I3OMGeoMap: Map of Optical Modules (OMs).
        
    Raises:
        Exception: If unknown optical module type is specified.
    """
    geomap = dataclasses.I3OMGeoMap()
    
    # Add OMs at x,y and depth
    for i in range(len(xpositions)):
        for m in range(len(depth)):
            loc = [xpositions[i], ypositions[i], depth[m]]
            string_num = i+1
            om_num = m+1
            if omsequence[m%len(omsequence)]=='POM':
                AddPOM(geomap, loc, string_num, om_num)
            elif omsequence[m%len(omsequence)]=='PCAL':
                AddPCAL(geomap, loc, string_num, om_num)
            else:
                raise Exception("Unknown optical module type specified. Type can be 'POM' or 'PCAL'")
    
    return geomap


def generateGFrame(xpositions, ypositions, depthlist, omsequence, domradius):
    """Creates a geometry frame with all the needed fields in a G frame.
    
    Args:
        xpositions (list): List of x positions for strings.
        ypositions (list): List of y positions for strings.
        depthlist (list): List of depths for OMs.
        omsequence (list): List defining the OM sequence (POM/PCAL).
        domradius (float): Radius of the optical modules.
        
    Returns:
        I3Frame: An I3Frame object to be used as the G frame.
    """
    geometry = dataclasses.I3Geometry()
    geomap = generateGeometry(xpositions, ypositions, depthlist, omsequence)
    geometry.omgeo = geomap

    # Create orientation - OM oriented upwards
    orientation = dataclasses.I3Orientation(0, 0, 1, 1, 0, 0)
    
    # Initialize frame as a G frame
    gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
    
    # Add geometry to frame
    gframe["I3Geometry"] = geometry
    #Duplicate copy of geomap in frame. Taking out for now
    #gframe["I3OMGeoMap"] = geomap
    
    # Create module geometry map
    modgeomap = dataclasses.I3ModuleGeoMap()
    for dom in geomap.keys():
        mkey = dataclasses.ModuleKey(dom.string, dom.om)
        module = dataclasses.I3ModuleGeo()
        module.module_type = dataclasses.I3ModuleGeo.ModuleType.mDOM
        module.orientation = orientation
        module.pos = dataclasses.I3Position(
            xpositions[dom.string-1], ypositions[dom.string-1], depthlist[dom.om-1]
        )
        module.radius = domradius
        modgeomap[mkey] = module
    gframe["I3ModuleGeoMap"] = modgeomap
    
    # Setting subdetector labels
    subdetec = dataclasses.I3MapModuleKeyString()
    for dom in geomap.keys():
        mkey = dataclasses.ModuleKey(dom.string, dom.om)
        subdetec[mkey] = omsequence[dom.om%len(omsequence)]
    gframe["Subdetectors"] = subdetec
    
    # Add start and end times
    gframe["StartTime"] = start_time
    gframe["EndTime"] = end_time
    
    return gframe
