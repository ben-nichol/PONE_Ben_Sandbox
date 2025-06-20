# written by jakub stacho based on code from thomas mcelroy
# last modified May 30, 2025

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

    om_times = {}
    om_pmts  = {}
    
    # extract from all individual pmt pulse maps
    for pmtkey in pulse_map.keys():
        # ignore this omkey if it is supposed to be dropped
        if pmtkey.string in drop_strings:
            continue
        if pmtkey.om in drop_oms:
            continue

        omkey = NoPMTKey(ModuleKey(pmtkey.string, pmtkey.om))
        if omkey not in om_times.keys():
            om_times[omkey] = []
        if omkey not in om_pmts.keys():
            om_pmts[omkey] = []

        for pulse in pulse_map[pmtkey]:
            om_times[omkey].append(pulse.time)
            om_pmts[omkey].append(pmtkey.pmt)

    for omkey in om_times.keys():
        sorted_om_mcpes = sorted(zip(om_times[omkey], om_pmts[omkey]))
        mcpe_map[omkey] = sorted_om_mcpes

    return mcpe_map