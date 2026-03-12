import os
import sys
import numpy as np
import scipy.constants as const
from os.path import dirname, abspath

from scipy.interpolate import CubicSpline, LinearNDInterpolator

from icecube import phys_services
from icecube.icetray import I3Units

class POM:
    '''
    Simple model of the POM
    '''
    def __init__(self,
                 qe_file        = os.getenv('PONESRCDIR') + '/data/qe.csv',
                 aa_file        = os.getenv('PONESRCDIR') + '/data/Angular_Acceptance.txt',
                 glass_file     = os.getenv('PONESRCDIR') + '/data/glass.csv',
                 random_service = None):
        '''
        Define P-OM properties and geometry
        '''
        # ----------------------------------------------------------------------------
        # set up pmt geometry
        # ----------------------------------------------------------------------------
        self.MODULE_RADIUS_M      = 0.2159
        # self.PMT_RADIUS_M         = 0.055
        self.PMT_RADIUS_M         = 0.1 # CHANGED THIS TO MATCH THE ANGULAR ACCEPTANCE FILES
                                        # This now allows 2 PMTs to be hit as the same time
                                        # Can choose what one is hit based on acceptance probability
        # [Zenith, Azimuth] angles for each PMT
        # indices in this array correspond to the
        # pmt number in the POM frame -1
        self.PMT_ANGLES = np.array([[  58.,0.  ],
                                    [  90., 328.  ],
                                    [ 122., 0.  ],
                                    [  90., 32.  ],
                                    [  51.37, 53.06],
                                    [  51.37, 306.94],
                                    [ 128.63, 306.94],
                                    [ 128.63, 53.06],
                                    [  58., 180.  ],
                                    [  90., 148.  ],
                                    [ 122., 180.  ],
                                    [  90., 212.  ],
                                    [  51.37, 233.06],
                                    [  51.37, 126.94],
                                    [ 128.63, 126.94],
                                    [ 128.63, 233.06]])

        # define xyz unit vectors for all PMT centres
        # PMTs 1-8 will be pointing in the +x direction and PMTs 9-16 will be pointing in the -x direction.
        x_coordinates = np.multiply(np.sin(np.deg2rad(self.PMT_ANGLES[:, 0])), np.cos(np.deg2rad(self.PMT_ANGLES[:, 1])))
        y_coordinates = np.multiply(np.sin(np.deg2rad(self.PMT_ANGLES[:, 0])), np.sin(np.deg2rad(self.PMT_ANGLES[:, 1])))
        z_coordinates = np.cos(np.deg2rad(self.PMT_ANGLES[:, 0]))

        self.PMT_MATRIX = np.array([x_coordinates, y_coordinates, z_coordinates]).T

        # ----------------------------------------------------------------------------
        # set up efficiencies
        # ----------------------------------------------------------------------------
        self.qe_function    = self.get_qe_function(qe_file)
        self.aa_function    = self.get_aa_function(aa_file)
        self.glass_function = self.get_glass_function(glass_file)

        self.collection_efficiency = 0.9 # assume a fixed constant collection efficiency for now

        # ----------------------------------------------------------------------------
        # set up random service if using icetray (otherwise just uses numpy random)
        # ----------------------------------------------------------------------------
        self.random_service = random_service

    # ----------------------------------------------------------------------------
    # methods for loading data
    # ----------------------------------------------------------------------------
    def get_qe_function(self, qe_file):
        '''
        Function that loads QE data
        '''
        energy_conversion = const.h * const.c / const.elementary_charge

        qe_data     = np.loadtxt(qe_file, delimiter=',')
        qe_function = CubicSpline(np.flip(energy_conversion / qe_data.T[0], axis=0), np.flip(qe_data.T[1], axis=0), extrapolate=False)

        return qe_function
    

    def get_aa_function(self, aa_file):
        '''
        Function that loads single PMT angular
        acceptance data
        '''
        aa_data     = np.loadtxt(aa_file, skiprows=1, delimiter=',')
        points = np.column_stack([aa_data.T[0], aa_data.T[1]])
        aa_function = LinearNDInterpolator(points, aa_data.T[2])

        return aa_function
    

    def get_glass_function(self, glass_file):
        '''
        Function that loads glass transmission
        data
        '''
        glass_data      = np.loadtxt(glass_file, delimiter=',')
        glass_function  = CubicSpline(glass_data.T[0] * 1e-9, glass_data.T[1], extrapolate=False)
    
        return glass_function


    # ----------------------------------------------------------------------------
    # methods for collecting efficiency values
    # ----------------------------------------------------------------------------
    def get_quantum_efficiency(self, photon_list):
        '''
        Returns an array of quantum efficiencies for
        a given array of photons
        '''
        wavelengths  = np.array([p.wavelength for p in photon_list])
        efficiencies = np.nan_to_num(self.qe_function(wavelengths))

        return efficiencies
    

    def get_angular_acceptance(self, pmt_distance_list, pmt_angle_list):
        '''
        Returns an array of angular acceptances
        for a given array of distances and angles between
        pmts and photons
        '''
        distances = np.atleast_1d(pmt_distance_list)
        angles    = np.atleast_1d(pmt_angle_list)
        points    = np.column_stack([distances, angles])  # shape (N, 2)
        
        efficiencies = np.nan_to_num(self.aa_function(points))
        return efficiencies

    def get_glass_efficiency(self, photon_list):
        '''
        Returns an array of glass transmittance efficiencies
        for a given array of photons
        '''
        wavelengths  = np.array([p.wavelength for p in photon_list])
        efficiencies = np.nan_to_num(self.glass_function(wavelengths))

        return efficiencies
    

    def get_collection_efficiency(self, photon_list):
        '''
        Returns an array of collection efficiencies
        for a given array of photons
        '''
        return np.ones_like(photon_list) * self.collection_efficiency


    # ----------------------------------------------------------------------------
    # methods for applying the full POM model
    # ----------------------------------------------------------------------------
    def get_pmt(self, photon):
        '''
        Returns the PMTs that are hit by a given
        photon based on geometry along with the
        distance and angle relative to the 
        corresponding PMT direction.
        '''
        photon_position            = photon.pos
        photon_position_normalized = np.array([photon_position.x,photon_position.y, photon_position.z]) / photon_position.magnitude
        photon_direction           = photon.dir
        photon_direction           = np.array([photon_direction.x, photon_direction.y, photon_direction.z])
         # deterine the angle between all the pmt vectors and the photon vector
        pmt_photon_angles = self.PMT_MATRIX.dot(photon_position_normalized)
        # check if the ditance at the module radius between the pmt vector
        # and the photon vector falls within the PMT size, if so, mark as a hit
        # return a list of indices that have been 'hit'
        # THIS SHOULD BE THE DISTANCE THAT IS ORTHOGONAL TO THE PMT VECTOR
        pmt_vector_distances = self.MODULE_RADIUS_M * np.sin(np.arccos(pmt_photon_angles))
        # Determine if the angle between the photon direction and the vector from the PMT to the photon is positive or negative 
        
        pmts_hit             = np.where(np.logical_and(pmt_vector_distances < self.PMT_RADIUS_M, pmt_photon_angles >= 0))[0]
        if len(pmts_hit) > 0:
            distance_vector = self.PMT_MATRIX[pmts_hit] - photon_position_normalized
            angle_sign = np.sign(distance_vector.dot(photon_direction))
            cos_angle = np.dot(self.PMT_MATRIX[pmts_hit], - photon_direction)
            # Will return multiple hits, need to apply acceptance cut to determine which one is actually detected
            hit_angle    = np.rad2deg(np.arccos(cos_angle)) * angle_sign
            hit_distance = pmt_vector_distances[pmts_hit]
            return np.array([[pmts_hit[i]+1, hit_distance[i], hit_angle[i]] for i in range(len(pmts_hit))])
        else:
            # no pmt hit return a PMT label 100 to mark an error
            return np.array([[100., 100., 100.]])


    def get_probabilities(self, photons, pmt_photon_distances,pmt_photon_angles):
        '''
        Returns a list of detection probabilities
        correspon to the given list of photons and
        distances between photons and PMTs
        '''
        weights = np.array([p.weight for p in photons])

        quantum_efficiency    = self.get_quantum_efficiency(photons)
        angular_acceptance    = self.get_angular_acceptance(pmt_photon_distances,pmt_photon_angles)
        glass_transmittance   = self.get_glass_efficiency(photons)
        collection_efficiency = self.get_collection_efficiency(photons)

        return weights * quantum_efficiency * angular_acceptance * glass_transmittance * collection_efficiency
    

    def module_acceptance(self, photon_list):
        '''
        Returns a list of photons, hit PMTs, and
        hit probabilities based on the module
        response to the given photon list
        '''
        pmt_list          = np.zeros_like(photon_list, dtype=float)
        hit_distance_list = np.zeros_like(photon_list, dtype=float)
        hit_angle_list    = np.zeros_like(photon_list, dtype=float)
        for i, photon in enumerate(photon_list): # NEED TO CHANGE THIS TO ACCOUNT FOR MULTIPLE PMT HITS
            
            hit_list = self.get_pmt(photon)
            
            angular_prob = self.get_angular_acceptance(hit_list.T[0], hit_list.T[1])
            best_hit = angular_prob==np.max(angular_prob)
            pmt_list[i]          = hit_list[best_hit][0][0]
            hit_distance_list[i] = hit_list[best_hit][0][1]
            hit_angle_list[i]    = hit_list[best_hit][0][2]

        probability_list = self.get_probabilities(photon_list, hit_distance_list,hit_angle_list)

        accepted_indices = np.logical_and(probability_list > 0.001, pmt_list != 100.)

        return [photon_list[accepted_indices], pmt_list[accepted_indices], probability_list[accepted_indices]]

    
    def apply_acceptance_cut(self, photon_list):
        '''
        Uses the module_acceptance function to
        attribute photons to PMT but also applies
        a cut based on the detection probability
        of each PMT. Returns a list of of surviving
        photons and PMT ids after the cut is applied
        '''
        photons, pmts, probabilities = self.module_acceptance(photon_list)
        if self.random_service is None:
            sampled_random_numbers = np.random.uniform(size = len(probabilities))
        else:
            sampled_random_numbers = np.array([self.random_service.uniform(0.0, 1.0) for i in np.arange(len(probabilities))])

        surviving_pmts    = pmts[sampled_random_numbers <= probabilities]
        surviving_photons = photons[sampled_random_numbers <= probabilities]

        return [surviving_photons, surviving_pmts]
    