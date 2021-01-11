#!/usr/bin/env python

'''
Pulse Cleaning method pulled from Dvir Hilu's original SimAnalysis.py script.

'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from numpy import linalg as la
from icecube.phys_services import I3Calculator

# A modified binary search algorithm to find all hits within a 20ns time 
# window. Takes entire list of MCPE hits, finds a 20ns window where the
# number of hits surpasses the hit threshold, and outputs those hits. If
# there are multiple windows, only the first occurence (in time) is taken. 
# NOTE: function assumes the list has already passed through the passDOM 
#       function. If the resulting list does not meet that criteria, the 
#       function raises a ValueError as it could cause issues for future 
#       functions relying on the significant MCPEs 
# 
# @Param:
# mcpeList:         The list of hits the DOM registered   
# hitThresh:        The number of hits required in the window for the DOM to
#                   be passed    
# 
# @Return: 
# An MCPESeries object with all the hits that were within the time window
def getSignificantMCPEs(mcpeList, hitThresh):
    timeList = [mcpe.time for mcpe in mcpeList]
    timeList.sort()
    highIndex = 0
    lowIndex = 0
    significantMCPEList = simclasses.I3MCPESeries()

    for mcpe in mcpeList:

        lowest = 0
        highest = len(timeList)-1
        while(lowest <= highest):
            middle = int((lowest+highest)/2)
            # if mcpe.time is greater than or equal to timeList[middle],  
            # then search in timeList[middle + 1 to highest]
            if(timeList[middle] <= mcpe.time + 20): 
                lowest = middle + 1 
            else: 
            # else search in timeList[llowest to middle-1] 
                highest = middle - 1
        
        highIndex = highest

        lowest = 0
        while(lowest <= highest):
            middle = int((lowest+highest)/2)
            # if mcpe.time is greater than or equal to timeList[middle],  
            # then search in timeList[middle + 1 to highest]
            if(timeList[middle] < mcpe.time): 
                lowest = middle + 1 
            else: 
            # else search in timeList[lowest to middle-1] 
                highest = middle - 1
        
        lowIndex = highest + 1

        if highIndex - lowIndex >= hitThresh - 1:
            break
    
    for mcpe in mcpeList:
        if mcpe.time >= timeList[lowIndex] and mcpe.time <= timeList[highIndex]:
            significantMCPEList.append(mcpe)

    # sanity check
    if len(significantMCPEList) < hitThresh:
        raise ValueError("There are not enough hits in this list. Make sure to filter with passDOM before calling this function")
    
    return significantMCPEList

# Converts the MCPESeriesMap to one that contains siginficant hits. Significant
# hits are defined to be ones where the DOM passed the passDOM function. The
# function then writes the new MCPESeriesMap to the frame before returning it
# NOTE: function assumes the list has already passed through the passFrame 
#       function. If the resulting map does not meet that criteria, the a 
#       ValueError is raised as it could cause issues for future functions 
#       relying on the significant MCPESeriesMap
#  
# @Param:
# frame:            The frame to be appended then returned   
# domsUsed:         The omkeys used for the analysis. Allows to look at
#                   smaller geometries from a larger geometry sim file
# hitThresh:        Hit threshold for a single DOM
# domThresh:        Number of DOMs passing the hit threshold to pass the
#                   frame  
# maxResidual:      Maximum time residual allowed for a hit to be considered
# GeoMap:           An I3OMGeoMap object that maps the omkeys in the
#                   geometries to their I3OMGeo objects 
# 
# @Return: 
# The frame inputted to the function after the significant MCPESeriesMap object
# was appended to it with the key word "MCPESerieMap_significant_hits" 
def writeSigHitsMapToFrame(frame, domsUsed, hitThresh, domThresh, geoMap):

    mcpeMap = frame["MCPESeriesMap"]
    significantMCPEMap = simclasses.I3MCPESeriesMap()
    for omkey in domsUsed:
        position = geoMap[omkey].position
        if omkey not in mcpeMap:
            continue
        if len(mcpeList) >= hitThresh :
            significantMCPEMap[omkey] = getSignificantMCPEs(mcpeMap[omkey], hitThresh)  
    
    # sanity check
    if len(mcpeMap) < domThresh:
        raise ValueError("Not enough DOMs in map. Make sure to filter with passFrame before calling this function")
    
    frame.Put("MCPESeriesMap_significant_hits", significantMCPEMap)

    return frame