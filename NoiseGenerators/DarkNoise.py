import numpy as np

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey

from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from NoiseGenerators.NoiseUtility import get_noise_time_bounds, get_mcpe_map

import sys
sys.path.append('/home/jakubs/projects/def-mdanning/jakubs/k40/utils')
from POMModel import POM



class DarkNoise(icetray.I3ConditionalModule):
    '''
    Icetray module that create a pulse map for
    dark noise
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('input_map',
                          'Name of the i3phtons after module acceptance is applied',
                          'I3Photons_pmtsplit')
        self.AddParameter('output_map',
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
        self.output_map         = self.GetParameter('output_map')
        self.dark_rate          = self.GetParameter('dark_rate')
        self.random_service     = self.GetParameter('random_service')
        self.drop_strings       = self.GetParameter('drop_strings')
        self.drop_oms           = self.GetParameter('drop_oms')
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
            mcpe_map  = get_mcpe_map(pulse_map, self.drop_strings, self.drop_oms)

            lower_time_limit, upper_time_limit = get_noise_time_bounds(mcpe_map, self.noise_padding[0], self.noise_padding[1])
        
        frame[self.output_map] = self.generate_dark_hits(lower_time_limit, upper_time_limit)

        self.PushFrame(frame)