from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
from Utilities.DOMUtility import NoPMTKey, AddPMTKey

import sys
sys.path.append('/home/jakubs/projects/def-mdanning/jakubs/k40/utils')
from POMModel import POM



class DarkNoise(icetray.I3ConditionalModule):
    '''
    Icetray module that separates I3Photons from CLSim into
    their respective PMTs and applies a POM acceptance cut.
    Writes the surviving photons back into the frame as I3RecoPulses
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('input_map',
                          'Name of the i3phtons after module acceptance is applied',
                          'I3Photons_pmtsplit')
        self.AddParameter('dark_name',
                          'Name for the tree of dark hits in the i3 file',
                          'DarkHits')
        self.AddParameter('dark_rate',
                          'Dark Noise rate (pulses per ns)',
                          0.000001)
        self.AddParameter('random_service',
                          'I3RandomService')
        self.AddParameter('drop_strings',
                          'List of string indices to ignore',
                          [])
        self.AddParameter('drop_oms',
                          'List of om indices to ignore',
                          [])
        self.AddParameter('drop_empty',
                          'Bool to determine if empty frames should be removed from the i3 file',
                          False)
        self.AddParameter('use_manual_noise_bounds',
                          'Bool to determine whether to use manual time bounds or calculate them per frame',
                          False)
        self.AddParameter('manual_noise_bounds',
                          'Manually specify the time bounds [min, max] over which to generate noise in ns',
                          [0., 10_000.])
        self.AddParameter('noise_padding',
                          'Specify the time padding bounds for noise [before, after] before the first and after the last photon in ns (Only applied if not using manual noise bounds)',
                          [2_000.0, 10_000.0])
        self.AddParameter('gcd_file',
                          'GCD file to use if these frames are not included in the i3 files',
                          None)


    def Configure(self):
        self.input_map          = self.GetParameter('input_map')
        self.dark_name          = self.GetParameter('dark_name')
        self.dark_rate          = self.GetParameter('dark_rate')
        self.random_service     = self.GetParameter('random_service')
        self.drop_strings       = self.GetParameter('drop_strings')
        self.drop_oms           = self.GetParameter('drop_oms')
        self.drop_empty         = self.GetParameter('drop_empty')
        self.use_manual_bounds  = self.GetParameter('use_manual_noise_bounds')
        self.manual_time_bounds = self.GetParameter('manual_noise_bounds')
        self.noise_padding      = self.GetParameter('noise_padding')

        gcd_file = self.GetParameter('gcd_file')
        if gcd_file is not None:
            # load in the geometry from a gcd file if required
            gcd_file = dataio.I3File(gcd_file, 'r')
            for frame in gcd_file:
                self.omkeys_to_use = frame['I3ModuleGeoMap'].keys()
                break
            gcd_file.close()
        else:
            self.omkeys_to_use = None
        
        self.num_pmts = 16
    

    def get_mcpe_map(self, pulse_map, drop_strings=[], drop_oms=[]):
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
    

    def get_noise_time_bounds(self, mcpe_map):
        '''
        Get the upper and lower time bounds to set a time
        window for dark hit generation.
        
        NOTE: If the mcpe map is empty, this will return a
        lower bound that is larger than the upper bound and
        so the generate_dark_hits method will return an empty
        pulse map
        '''
        lower_bound = 999999999.0
        upper_bound = -999999999.0

        for omkey in mcpe_map.keys():
            for pulse in mcpe_map[omkey]: # pulse is a tuple (time, pmt)
                if upper_bound < pulse[0]:
                    upper_bound = pulse[0]
                if lower_bound > pulse[0]:
                    lower_bound = pulse[0]

        lower_bound -= self.noise_padding[0]
        upper_bound += self.noise_padding[1]

        return lower_bound, upper_bound


    def generate_dark_hits(self, lower_bound, upper_bound):
        '''
        Add dark hits across all OMs (Adds I3MCPEs)
        '''
        dark_pulse_map = dataclasses.I3RecoPulseSeriesMap()

        if lower_bound > upper_bound:
            return dark_pulse_map
        
        for omkey in self.omkeys_to_use:
            if omkey.string in self.drop_strings:
                continue
            if omkey.om in self.drop_oms:
                continue

            for pmt_index in np.arange(self.num_pmts):
                pmtkey = AddPMTKey(omkey, pmt_index)

                # poisson distribution of dark noise so sample
                # time between hits from an exponential distribution
                time_delta = self.random_service.exp(1.0 / self.dark_rate) + lower_bound

                first = True
                while time_delta < upper_bound:
                    if first:
                        dark_pulse_map[pmtkey] = dataclasses.I3RecoPulseSeries()
                        first = False

                    pulse        = dataclasses.I3RecoPulse()
                    pulse.time   = time_delta
                    pulse.charge = 1.0
                    dark_pulse_map[pmtkey].append(pulse)

                    # get the time to the next dark hit
                    # and the next PMT to be hit
                    time_delta += self.random_service.exp(1.0 / self.dark_rate)
        
        return dark_pulse_map


    def Geometry(self, frame):
        if self.omkeys_to_use is None:
            self.omkeys_to_use = frame['I3Geometry'].omgeo.keys()

        self.PushFrame(frame)


    def DAQ(self, frame):
        if self.use_manual_bounds:
            lower_time_limit = self.manual_time_bounds[0]
            upper_time_limit = self.manual_time_bounds[1]
        else:
            pulse_map = frame[self.input_map]
            mcpe_map  = self.get_mcpe_map(pulse_map, self.drop_strings, self.drop_oms)

            lower_time_limit, upper_time_limit = self.get_noise_time_bounds(mcpe_map)
        
        frame[self.dark_name] = self.generate_dark_hits(lower_time_limit, upper_time_limit)

        self.PushFrame(frame)