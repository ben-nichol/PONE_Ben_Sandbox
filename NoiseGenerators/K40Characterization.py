class K40Characterization:
    '''
    Class that stores the characterization
    information so we can save as a pickle

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
    def __init__(self, combinations, probabilities, multifold_rate, single_rate, arrival_time_spread):
        self.combinations        = combinations
        self.probabilities       = probabilities
        self.multifold_rate      = multifold_rate
        self.single_rate         = single_rate
        self.arrival_time_spread = arrival_time_spread
        

        self.single_fraction = single_rate / (multifold_rate + single_rate)
        self.multi_fraction  = 1 - self.single_fraction
        self.total_rate_ns   = (multifold_rate + single_rate) / 1e9