#!/usr/bin/env python

"""
a few helper functions that will be commonly used when making
GCD files to test prototype geometries
"""

from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units
from enum import Enum
import numpy as np

cdfile = dataio.I3File("Calib_and_DetStat_File.i3.gz")
cdframe = cdfile.pop_frame()
calib = cdframe["I3Calibration"]
start_time = calib.start_time
end_time = calib.end_time
bedrockz = -I3Constants.OriginElev


# Takes the distance from the surface and returns the z position in
# Icecube coordinates.
#
# @Param:
# depth:        a float representing vertical distance from the surface
#
# @Return:
# the corresponding z value in the IceCube Coordinate system
def convertDepthToZ(depth):
    return I3Constants.SurfaceElev - I3Constants.OriginElev - depth


# Creates a calibration with all the correct OMKeys
#
# @Param:
# geometry:     the I3Geometry object inputted to the frame
#
# @Return:
# An I3Calibration object to be used in the GCD file
def makeCalibrationObject(geometry):
    calibrationData = cdframe["I3Calibration"].dom_cal.values()[0]
    domcal = dataclasses.Map_OMKey_I3DOMCalibration()
    calibration = cdframe["I3Calibration"]
    for omkey in geometry.omgeo.keys():
        domcal[omkey] = calibrationData

    calibration.dom_cal = domcal

    return calibration


# Creates a detector status with all the correct OMKeys
#
# @Param:
# geometry:     the I3Geometry object inputted to the frame
#
# @Return:
# An I3DetectorStatus object to be used in the GCD file
def makeDSObject(geometry):
    dsData = cdframe["I3DetectorStatus"].dom_status.values()[0]
    domstat = dataclasses.Map_OMKey_I3DOMStatus()
    detectorStatus = cdframe["I3DetectorStatus"]
    for omkey in geometry.omgeo.keys():
        domstat[omkey] = dsData

    detectorStatus.dom_status = domstat

    return detectorStatus


# Creates a frame with all the needed fields in a C frame
#
# @Param:
# geometry:     the I3Geometry object inputted to the frame
#
# @Return:
# An I3Frame object to be used as the C frame
def generateCFrame(geometry):
    # intialize frame as a C frame
    frame = icetray.I3Frame(icetray.I3Frame.Calibration)

    # add key-value pairs
    frame["I3Geometry"] = geometry
    frame["I3Calibration"] = makeCalibrationObject(geometry)
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]

    return frame


# Creates a frame with all the needed fields in a D frame
#
# @Param:
# geometry:     the I3Geometry object inputted to the frame
#
# @Return:
# An I3Frame object to be used as the D frame
def generateDFrame(geometry):
    # intialize frame as a D frame
    frame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

    # add key-value pairs
    frame["I3Geometry"] = geometry
    frame["I3Calibration"] = makeCalibrationObject(geometry)
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]
    frame["I3DetectorStatus"] = makeDSObject(geometry)
    frame["BadDomsList"] = cdframe["BadDomsList"]
    frame["BadDomsListSLC"] = cdframe["BadDomsListSLC"]

    return frame


# Generates a string of DOMs in a specified direction.
# This is represented by an I3OMGeoMap object, which
# has an OMKey object for a DOM and an I3OMGeo object
# to represent its geometry.
#
# @Param:
# stringNumber: an integer representing the string id
# startPos:     an I3Position object with the location
#               of the top of the string
# numDoms:      number of DOMs on the string
# spacing:      array detailing DOM spacing on string
# direction:    an I3Direction object representing the
#               direction of the line of the string
#
# @Return:
# an I3OMGeoMap object with the DOMs on the string and their
# geometries
def generateOMString(stringNumber, startPos, numDoms, spacing, direction):
    orientation = dataclasses.I3Orientation(
        0, 0, -1, 1, 0, 0
    )  # same orientation as icecube DOMs (dir=down)
    area = 0.5857538 * I3Units.meter2  # same area as KM3NET MDOMs
    geomap = dataclasses.I3OMGeoMap()
    x = startPos.x
    y = startPos.y
    z = startPos.z
    dx = [spacingVal * direction.x for spacingVal in spacing]
    dy = [spacingVal * direction.y for spacingVal in spacing]
    dz = [spacingVal * direction.z for spacingVal in spacing]

    # create OMKeys and I3OMGeo for DOMs on string and add them to the map
    for i in xrange(0, numDoms):
        omkey = OMKey(stringNumber, i, 0)
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.IceCube
        omGeometry.orientation = orientation
        omGeometry.area = area
        omGeometry.position = dataclasses.I3Position(
            x + dx[i] * i, y + dy[i] * i, z + dz[i] * i
        )
        geomap[omkey] = omGeometry

    return geomap


# Generates a line of DOMs with vertical strings. This was
# necessary due to complications caused by clsim's xy divisions.
# The simulation divides the grid into x-y cells, but this
# is only compatible with verticle strings since it tries to
# divide it such that only one string fits in an x-y cell.
#
# @Param:
# stringNum:    an integer representing the id of the first
#               string in the DOM line
# startPos:     an I3Position object with the location of
#               the start of the DOM line
# spacing:      array detailing DOM spacing on line
# direction:    an I3Direction object representing the
#               direction of the line
# vertSpacing:  integer representing the distance between layers
# numStrings:   number of strings per DOM line
# layers:       number of DOMs in each string in the DOM line
#
# @Return:
# an I3GeoMap objects with the geometries of DOMs in the line
def generateDOMLine(
    stringNum, startPos, spacing, direction, vertSpacing, numStrings, layers
):
    lineMap = dataclasses.I3OMGeoMap()
    x = startPos.x
    y = startPos.y
    z = startPos.z
    dx = [spacingVal * direction.x for spacingVal in spacing]
    dy = [spacingVal * direction.y for spacingVal in spacing]
    dz = [spacingVal * direction.z for spacingVal in spacing]
    stringSpacing = [vertSpacing for i in range(0, layers)]

    for i in xrange(0, numStrings):
        currentNum = stringNum + i
        stringPos = dataclasses.I3Position(x + dx[i] * i, y + dy[i] * i, z + dz[i] * i)
        stringDirection = dataclasses.I3Direction(0, 0, 1)
        stringMap = generateOMString(
            currentNum, stringPos, layers, stringSpacing, stringDirection
        )
        lineMap.update(stringMap)

    return lineMap


# Generates a list of offsets to the DOM string starting positions. Different
# offset types are described in OffsetType enum class
#
# @Param:
# offset_type:  an enum indicating which offset method is chosen
# length:       the length of the resulting list. Should be equal to number
#               of strings in the layer
#
# @Return:
# a list containing different offset values for the starting positions of the
# DOM strings
def generateOffsetList(offset_type, length):
    offsetList = []

    if not isinstance(offset_type, OffsetType):
        raise TypeError("offset_type must be an instance of OffsetType Enum")

    if offset_type == OffsetType.LinearResetOffset:
        offset = 0
        for i in xrange(0, length):
            if offset > 300:
                offset = 100
            offsetList.append(offset)
            offset += 100
    elif offset_type == OffsetType.LinearRiseFallOffset:
        # start at center
        offset = 0
        # determines whether rising or falling offset
        signFactor = 1
        for i in xrange(0, length):
            if offset >= 300:
                signFactor = -1
            if offset <= 150:
                signFactor = 1
            offsetList.append(offset)
            offset += 150 * signFactor
    else:
        offsetList = [50 for i in range(0, length)]

    return offsetList


# Generates a list detailing the spacings between DOMs along a string. Different
# spacing types are described in SpacingType enum class
#
# @Param:
# spacing_type:     an enum indicating which spacing method was chosen
# basicSpacing:     the spacing between the first two DOMs
# length:           the length of the resulting list. Should be equal to number
#                   of DOMs in the string
#
# @Return:
# a list containing different spacing values for the DOMs along the string
def generateSpacingList(spacing_type, basicSpacing, length):
    spacingList = []
    undistortedStringLength = basicSpacing * length

    if not isinstance(spacing_type, SpacingType):
        raise TypeError("spacing_type must be an instance of SpacingType Enum")

    if spacing_type == SpacingType.LinearRSpacing:
        r = 0
        for i in xrange(0, length):
            spacing = basicSpacing * (1 - (r / undistortedStringLength))
            spacingList.append(spacing)
            r += spacing
    elif spacing_type == SpacingType.ExpRSpacing:
        r = 0
        for i in xrange(0, length):
            spacing = basicSpacing * np.exp(-r / undistortedStringLength)
            spacingList.append(spacing)
            r += spacing
    else:
        spacingList = [basicSpacing for i in range(0, length)]

    print(sum(spacingList))
    return spacingList


# decomposes a geometry into a new geometry with only the DOMs specified
#
# @Param:
# denseGeometry:    an I3Geometry object representing the dense geometry
# domsUsed:         a list of OMKey objects representing the DOMs to be
#                   included in the new geometry. If an OMKey in domsUsed
#                   is not in denseGeometry, a ValueError is raised
#
# @Return:
# a new I3Geometry object containing only the DOMs in domsUsed
def makePartialGeometry(denseGeometry, domsUsed):
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


# an enum class to keep track of different offset types
#
# SimpleOffset:         offset of 50m to avoid DOMs overlapping
# LinearResetOffset:    Offset starts at and increases by 20m with each string.
#                       If it reaches 100m, it resets back to 20m.
# LinearRiseFallOffset: Offset starts at 0 and increases by 20m with each string.
#                       If it reaches 100m, offset starts decreasing by 20m.
#                       When reaching 20m, starts increasing again.
class OffsetType(Enum):
    SimpleOffset = "simple_offset"
    LinearResetOffset = "linear_reset_offset"
    LinearRiseFallOffset = "rise_fall_offset"

    def __str__(self):
        return self.value


# an enum class to keep track of different spacing types
#
# SimpleSpacing:        uniform spacing determined by the basicSpacing parameter
# LinearRSpacing:       Spacing decreases linearly with r according to the equation
#                       spacing = basicSpacing * ( 1 - (r/totalStringLength) )
# ExpRSpacing:          Spacing decreases exponentially with r according to the euqatiion
#                       spacing = basicSpacing * np.exp(-r/undistortedStringLength)
class SpacingType(Enum):
    SimpleSpacing = "simple_spacing"
    LinearRSpacing = "linear_r_spacing"
    ExpRSpacing = "exp_r_spacing"

    def __str__(self):
        return self.value
