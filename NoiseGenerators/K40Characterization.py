# written by jakub stacho
# last modified May 30, 2025

import h5py
import numpy as np

class K40Characterization:
    '''
    Class that stores the characterization information so we can
    access it easily. This class loads all the data from the
    k40 characterization hdf5 file

    combinations        - array of pmt coincidence combinations
    probabilities       - array of corresponding probabilities
    multifold_rate      - rate (Hz) of multi-fold coincidences
    single_rate         - rate (Hz) of single-fold coincidences
    arrival_time_spread - array of the difference in photon
                          arrival times(ns) [column 0] and the
                          corresponding probability [column 1]

    single_fraction     - fraction of coincidences that are single-fold
    multi_fraction      - fraction of coincidences that are multi-fold
    total_rate_ns       - the rate of all coincidences in 1/ns
    '''
    def __init__(self, characterization_file_path):
        with h5py.File(characterization_file_path, 'r') as h5f:
            # revert combinations back to a numpy
            # array of different sized numpy arrays
            padded_combinations      = h5f['coincidence-combinations/combinations'][:]
            self.combinations        = np.array([np.array(combination[combination >= 0]) for combination in padded_combinations], dtype=object)
            self.probabilities       = h5f['coincidence-combinations/weights'][:]
            self.single_rate         = h5f['coincidence-combinations'].attrs['singlefold-rate']
            self.multifold_rate      = h5f['coincidence-combinations'].attrs['multifold-rate']
            self.arrival_time_spread = h5f['arrival-times'][:]
        
            self.single_fraction = self.single_rate / (self.multifold_rate + self.single_rate)
            self.multi_fraction  = 1 - self.single_fraction
            self.total_rate_ns   = (self.multifold_rate + self.single_rate) / 1e9