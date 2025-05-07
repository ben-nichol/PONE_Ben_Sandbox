import sys
import numpy as np
import _pickle as pickle

from K40Characteriztaion import K40Characteriztaion

sys.path.append('/home/jakubs/projects/def-mdanning/jakubs/k40/utils/')
from POMModel import POM

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey



class Salt():
    '''
    Class for the P-ONE 40K background noise generator
    '''
    def __init__(self, characterization_file):
        swag
    

    def get_decay_times(self, frame_length_ns, event_rate):
        '''
        Returns a list of decay times within a frame of a
        given length using the expected event rate in (events/ns)
        '''
        mean    = frame_length_ns * event_rate
        samples = int(mean * 5)

        time_offsets    = np.random.gamma(1, 1/event_rate, samples)
        cumulative_time = np.cumsum(time_offsets)

        while True:
            try:
                # if the cumulative time has exceeded the frame length
                last_event_index = np.where(cumulative_time > frame_length_ns)[0][0]
                return cumulative_time[:last_event_index]
            except:
                # otherwise need to generate more times
                additional_offsets = np.random.gamma(1, 1/event_rate, samples)
                time_offsets       = np.concatenate((time_offsets, additional_offsets))
                cumulative_time    = np.cumsum(time_offsets)
    

    def get_photon_times(self, initial_time, num_photons, photon_time_distribution):
        '''
        Returns an array of photon arrival times based
        on the number of photons and initial decay time
        '''
        offsets = np.random.choice(photon_time_distribution[:,0], num_photons, p=photon_time_distribution[:,1])
        offsets += initial_time + np.random.uniform(-0.5, 0.5, num_photons)
        return offsets
    

    def distribute_pmts(self, generated_pmts, module, hflip, vflip, rotation):
        '''
        Function that takes a list of generated PMT indices
        including a home PMT and applies a random rotation
        and flip to them so that they can be evenly distributed
        over all the PMTs on the module

        hflip   : boolean for horizontal flip
        vflip   : boolean for vertical flip
        rotation: int 0-3 representing 0, 90, 180, 270 degree rotation
        '''
        generated_pmt_angles = module.PMT_ANGLES[generated_pmts]

        if module.UPPER_HOME_INDEX in generated_pmts:
            home = module.UPPER_RING_HOME
        elif module.LOWER_HOME_INDEX in generated_pmts:
            home = module.LOWER_RING_HOME
        else:
            raise Exception('No home PMT found in PMT list!')
        
        # horizontal flip
        if hflip:
            azimuth_differences = generated_pmt_angles[:, 1] - home[1]
            flipped_azimuths    = generated_pmt_angles[:, 1] - (2 * azimuth_differences)

            generated_pmt_angles[:, 1] = flipped_azimuths
            generated_pmt_angles[:,1][np.where(generated_pmt_angles[:, 1] < 0)[0]] += 360 # for plotting

        # vertical flip
        if vflip:
            generated_pmt_angles[:,0] *= -1
        
        # rotation
        rotation_angle = 90 * rotation
        generated_pmt_angles[:,1] = generated_pmt_angles[:,1] - rotation_angle
        generated_pmt_angles[:,1][np.where(generated_pmt_angles[:, 1] < 0)[0]] += 360

        distributed_pmts       = np.zeros(len(generated_pmts), dtype=int)
        distributed_angle_sums = np.sum(generated_pmt_angles, axis=1)

        for i, angle_sum in enumerate(distributed_angle_sums):
            distributed_pmts[i] = np.where(module.ANGLE_SUMS == angle_sum)[0]
        
        return distributed_pmts
    

    def sample_pmts(self, characterization, num_to_generate, photon_times, module):
        '''
        Samples pmt combinations and distributes pmts across the module.
        Also calculates the hit time at each PMT
        '''
        sampled_pmt_combinations = np.random.choice(characterization.combinations, num_to_generate, p=characterization.probabilities)

        distributed_pmts = np.empty(num_to_generate, dtype=object)
        pmt_hit_times    = np.empty(num_to_generate, dtype=object)

        hflips  = np.random.randint(0, 2, num_to_generate)
        vflips  = np.random.randint(0, 2, num_to_generate)
        rotates = np.random.randint(0, 4, num_to_generate)

        for i, pmts in enumerate(sampled_pmt_combinations):
            distributed_pmts[i] = self.distribute_pmts(pmts, module, hflips[i], vflips[i], rotates[i])
            pmt_hit_times[i]    = self.get_photon_times(photon_times[i], len(pmts), characterization.arrival_time_spread)
        
        return (np.hstack(distributed_pmts), np.hstack(pmt_hit_times))
    

    def get_k40(self, mkey, pulse_series_map, characterization, frame_time_ns, module):
        '''
        Main function for generating the k40
        and populating the pulse series
        '''
        coincidence_times = self.get_decay_times(frame_time_ns, characterization.total_rate_ns)

        # determine which coincidences are single fold and which are multi fold
        single_indices = np.random.choice([0, 1], len(coincidence_times), p=[characterization.multi_fraction, characterization.single_fraction]).astype(bool)
        multi_indices  = np.invert(single_indices)

        num_single   = np.sum(single_indices)
        single_times = coincidence_times[single_indices]
        single_pmts  = np.random.randint(0, 16, num_single)

        num_multi = np.sum(multi_indices)

        if num_multi == 0:
            pmts  = single_pmts
            times = single_times
        else:
            multi_pmts, multi_times = self.sample_pmts(characterization, num_multi, coincidence_times[multi_indices], module)

            pmts  = np.hstack((single_pmts, multi_pmts))
            times = np.hstack((single_times, multi_times))
        
        # remove any hits that fall outside of the frame time
        # window. This could occur because of the photon arrival
        # time spread that is applied in SamplePMTs
        out_of_bounds_time_indices = np.where(times > frame_time_ns)[0]
        if len(out_of_bounds_time_indices) > 0:
            times = np.delete(times, out_of_bounds_time_indices)
            pmts  = np.delete(pmts, out_of_bounds_time_indices)

        # conver from k40 numbering to poneoffline numbering
        pmts = module.k40_to_offline_pmts[pmts]

        string         = mkey.string
        optical_module = mkey.om

        # only add the omkeys for pmts that were hit to the pulse series
        for pmt in np.unique(pmts):
            pulse_series_map[OMKey(string, optical_module, int(pmt))] = dataclasses.I3RecoPulseSeries()

        for i, pmt in enumerate(pmts):
            pulse        = dataclasses.I3RecoPulse()
            pulse.time   = times[i]
            pulse.charge = 1.0

            pulse_series_map[OMKey(string, optical_module, int(pmt))].append(pulse)
            
        return pulse_series_map