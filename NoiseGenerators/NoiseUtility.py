from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import NoPMTKey, AddPMTKey



def get_noise_time_bounds(mcpe_map, noise_time_padding_before, noise_time_padding_after):
    '''
    Get the upper and lower time bounds to set a time
    window for noise generation.
    
    NOTE: If the mcpe map is empty, this will return a
    lower bound that is larger than the upper bound which
    should eventually result in an empty pulse map
    '''
    lower_bound = 999999999.0
    upper_bound = -999999999.0

    for omkey in mcpe_map.keys():
        for pulse in mcpe_map[omkey]: # pulse is a tuple (time, pmt)
            if upper_bound < pulse[0]:
                upper_bound = pulse[0]
            if lower_bound > pulse[0]:
                lower_bound = pulse[0]

    lower_bound -= noise_time_padding_before
    upper_bound += noise_time_padding_after

    return lower_bound, upper_bound



def get_mcpe_map(pulse_map, drop_strings=[], drop_oms=[]):
    '''
    Read a split pmt pulse map from the frame and
    return an OM wide mcpe map of hit times and PMTs
    '''
    mcpe_map = {}
    
    # make new map with individual PMTs
    for pmtkey in pulse_map.keys():
        # ignore this omkey if it is supposed to be dropped
        if pmtkey.string in drop_strings:
            continue
        if pmtkey.om in drop_oms:
            continue

        omkey = NoPMTKey(ModuleKey(pmtkey.string, pmtkey.om))
        if omkey not in mcpe_map.keys():
            mcpe_map[omkey] = []

        for pulse in pulse_map[pmtkey]:
            mcpe_map[omkey].append((pulse.time, pmtkey.pmt)) # mcpe map entries are tuples (time, pmt)

    return mcpe_map