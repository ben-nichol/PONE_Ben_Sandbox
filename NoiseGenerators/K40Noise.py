# written by jakub stacho
# last modified May 30, 2025

import os
import numpy as np

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey

from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.POMModel import POM
from NoiseGenerators.NoiseUtility import get_noise_time_bounds, get_mcpe_map
from NoiseGenerators.K40Characterization import K40Characterization



class K40Noise(icetray.I3ConditionalModule):
    '''
    Icetray module that create a pulse map for
    k40 noise
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('input_map',
                          'Name of the i3phtons after module acceptance is applied',
                          'I3Photons_pmtsplit')
        self.AddParameter('output_map',
                          'Name for the tree of k40 hits in the i3 file',
                          'K40Hits')
        self.AddParameter('characterization_file',
                          'Path to the k40 characterization hdf5 file',
                          os.getenv('PONESRCDIR') + '/NoiseGenerators/k40-characterization.hdf5')
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
        self.input_map             = self.GetParameter('input_map')
        self.output_map            = self.GetParameter('output_map')
        self.characterization_file = self.GetParameter('characterization_file')
        self.random_service        = self.GetParameter('random_service')
        self.drop_strings          = self.GetParameter('drop_strings')
        self.drop_oms              = self.GetParameter('drop_oms')
        self.use_manual_bounds     = self.GetParameter('use_manual_noise_bounds')
        self.manual_time_bounds    = self.GetParameter('manual_noise_bounds')
        self.noise_padding         = self.GetParameter('noise_padding')

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

        # load the updated module acceptance
        self.module   = POM()
        self.num_pmts = len(self.module.PMT_MATRIX)

        self.DEBUG_FRAME = 0


    def random_choice(self, values, probabilities, number_samples):
        '''
        Reproduces numpy.random.choice to sample from a
        given arbitrary probability histogram but using
        an I3RandomService rather than numpy random.
        Assumes probabilities are passed as a numpy array
        '''
        # make sure probabilities are normalized
        probabilities  /= np.sum(probabilities)
        probability_cdf = np.cumsum(probabilities)

        sampled_values = np.empty(number_samples, dtype=type(values[0]))
        for i in range(number_samples):
            sampled_index     = np.searchsorted(probability_cdf, self.random_service.uniform(0.0, 1.0))
            sampled_values[i] = values[sampled_index]
        
        return sampled_values


    def get_event_times(self, event_rate, lower_bound, upper_bound):
        '''
        Returns an array of coincidence times within
        a frame given the expected event rate (events/ns)
        and lower and upper time bounds (ns)
        '''
        event_times = []
        time_delta  = self.random_service.exp(1.0 / event_rate) + lower_bound

        while time_delta < upper_bound:
            event_times.append(time_delta)

            time_delta += self.random_service.exp(1.0 / event_rate)
        
        return np.array(event_times)
    

    def get_photon_times(self, initial_time, num_photons, photon_time_distribution):
        '''
        Returns an array of photon arrival times based on the
        number of photons and initial event time
        '''
        time_offset_bins          = photon_time_distribution[:, 0]
        time_offset_probabilities = photon_time_distribution[:, 1]
        half_time_offst_bin_size  = (time_offset_bins[1] - time_offset_bins[0]) / 2.

        time_offsets = self.random_choice(time_offset_bins, time_offset_probabilities, num_photons)

        # smear uniformly within the time bins
        for i in range(num_photons):
            bin_smear        = self.random_service.uniform(-half_time_offst_bin_size, half_time_offst_bin_size)
            time_offsets[i] += bin_smear

        return time_offsets + initial_time
    

    def distribute_pmts(self, generated_pmts, hflip, vflip, rotation):
        '''
        Function that takes a list of generated PMT indices
        including a home PMT and applies a random rotation
        and flip to them so that they can be evenly distributed
        over all the PMTs on the module

        hflip   : boolean for horizontal flip
        vflip   : boolean for vertical flip
        rotation: int 0-3 representing 0, 90, 180, 270 degree rotation
        '''
        generated_pmt_angles = self.module.PMT_ANGLES[generated_pmts]

        if self.module.UPPER_HOME_INDEX in generated_pmts:
            home = self.module.UPPER_RING_HOME
        elif self.module.LOWER_HOME_INDEX in generated_pmts:
            home = self.module.LOWER_RING_HOME
        else:
            raise Exception('No home PMT found in PMT list!')
        
        # horizontal flip
        if hflip:
            azimuth_differences = generated_pmt_angles[:, 1] - home[1]
            flipped_azimuths    = generated_pmt_angles[:, 1] - (2 * azimuth_differences)

            generated_pmt_angles[:, 1] = flipped_azimuths
            # make sure the angles are within [0, 360) degrees
            # so we can match with the defined pmt angles
            generated_pmt_angles[:, 1][np.where(generated_pmt_angles[:, 1] < 0)[0]] += 360
            generated_pmt_angles[:, 1][np.where(generated_pmt_angles[:, 1] >= 360)[0]] -= 360

        # vertical flip
        if vflip:
            generated_pmt_angles[:,0] *= -1
        
        # rotation
        rotation_angle = 90 * rotation
        generated_pmt_angles[:,1] = generated_pmt_angles[:,1] - rotation_angle
        
        # make sure the angles are within [0, 360) degrees
        # so we can match with the defined pmt angles
        generated_pmt_angles[:, 1][np.where(generated_pmt_angles[:, 1] < 0)[0]] += 360
        generated_pmt_angles[:, 1][np.where(generated_pmt_angles[:, 1] >= 360)[0]] -= 360

        distributed_pmts       = np.zeros(len(generated_pmts), dtype=int)
        distributed_angle_sums = np.sum(generated_pmt_angles, axis=1)

        for i, angle_sum in enumerate(distributed_angle_sums):
            distributed_pmts[i] = np.where(self.module.ANGLE_SUMS == angle_sum)[0][0]
        
        return distributed_pmts
    

    def sample_pmts(self, characterization, num_to_generate, photon_times):
        '''
        Samples pmt combinations and distributes pmts across the module.
        Also calculates the hit time at each PMT
        '''
        sampled_pmt_combinations = self.random_choice(characterization.combinations, characterization.probabilities, num_to_generate)

        distributed_pmts = np.empty(num_to_generate, dtype=object)
        pmt_hit_times    = np.empty(num_to_generate, dtype=object)

        hflips  = np.array([self.random_service.integer(2) for i in range(num_to_generate)])
        vflips  = np.array([self.random_service.integer(2) for i in range(num_to_generate)])
        rotates = np.array([self.random_service.integer(4) for i in range(num_to_generate)])

        for i, pmts in enumerate(sampled_pmt_combinations):
            distributed_pmts[i] = self.distribute_pmts(pmts, hflips[i], vflips[i], rotates[i])
            pmt_hit_times[i]    = self.get_photon_times(photon_times[i], len(pmts), characterization.arrival_time_spread)
        
        return (np.hstack(distributed_pmts), np.hstack(pmt_hit_times))
    

    def generate_k40_hits(self, characterization, lower_bound, upper_bound):
        '''
        Generate the k40 noise pulse map
        '''
        # Added changes to return I3MCPESeriesMap
        k40_mcpe_map = simclasses.I3MCPESeriesMap()

        if lower_bound > upper_bound:
            return k40_mcpe_map

        for omkey in self.omkeys_to_use:
            if omkey.string in self.drop_strings:
                continue
            if omkey.om in self.drop_oms:
                continue

            coincidence_times = self.get_event_times(characterization.total_rate_ns, lower_bound, upper_bound)

            # determine which coincidences are single fold and which are multi fold
            single_indices = self.random_choice([0, 1], [characterization.multi_fraction, characterization.single_fraction], len(coincidence_times)).astype(bool)
            multi_indices  = np.invert(single_indices)

            num_single   = np.sum(single_indices)
            single_times = coincidence_times[single_indices]
            single_pmts  = np.array([self.random_service.integer(self.num_pmts) for i in range(num_single)])

            num_multi = np.sum(multi_indices)

            if num_multi == 0:
                pmts  = single_pmts
                times = single_times
            else:
                multi_pmts, multi_times = self.sample_pmts(characterization, num_multi, coincidence_times[multi_indices])

                pmts  = np.hstack((single_pmts, multi_pmts))
                times = np.hstack((single_times, multi_times))
            
            # remove any hits that fall outside of the frame time
            # window. This could occur because of the photon arrival
            # time spread that is applied in sample_pmts
            out_of_bounds_time_indices = np.where(times > upper_bound)[0]
            if len(out_of_bounds_time_indices) > 0:
                times = np.delete(times, out_of_bounds_time_indices)
                pmts  = np.delete(pmts, out_of_bounds_time_indices)

            # only add the omkeys for pmts that were hit to the pulse series
            for pmt in np.unique(pmts):
                k40_mcpe_map[OMKey(omkey.string, omkey.om, int(pmt)+1)] = dataclasses.I3RecoPulseSeries()

            for i, pmt in enumerate(pmts):
                pulse        = simclasses.I3MCPE()
                pulse.time   = times[i]
                pulse.npe = 1 # Number of MCPEs
                # pulse.width  = int(pmt) # width seems to be a proxy for PMT number in the current DOMTrigger

                k40_mcpe_map[OMKey(omkey.string, omkey.om, int(pmt)+1)].append(pulse)

        return k40_mcpe_map


    def Geometry(self, frame):
        if self.omkeys_to_use is None:
            self.omkeys_to_use = frame['I3Geometry'].omgeo.keys()

        self.PushFrame(frame)


    def DAQ(self, frame):
        # print(f'FRAME # {self.DEBUG_FRAME}')
        # self.DEBUG_FRAME += 1
        # load k40 characterization
        characterization = K40Characterization(self.characterization_file)

        # determine noise generation time bounds
        if self.use_manual_bounds:
            lower_time_limit = self.manual_time_bounds[0]
            upper_time_limit = self.manual_time_bounds[1]
        else:
            pulse_map = frame[self.input_map]
            mcpe_map  = get_mcpe_map(pulse_map, self.drop_strings, self.drop_oms)

            lower_time_limit, upper_time_limit = get_noise_time_bounds(mcpe_map, self.noise_padding[0], self.noise_padding[1])
        
        frame[self.output_map] = self.generate_k40_hits(characterization, lower_time_limit, upper_time_limit)

        self.PushFrame(frame)
