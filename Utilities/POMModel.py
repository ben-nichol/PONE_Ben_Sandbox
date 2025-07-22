'''
Class defining an acceptance model for the POM
used for the K40 analysis and maybe other simulation
in the future
'''

import os
import sys
import numpy as np
import scipy.constants as const
from os.path import dirname, abspath

from scipy.interpolate import CubicSpline, interp1d

from icecube import phys_services
from icecube.icetray import I3Units



class POM:
    '''
    Simple model of the POM
    '''
    def __init__(self,
                 qe_file        = os.getenv('PONESRCDIR') + '/data/qe.csv',
                 aa_file        = os.getenv('PONESRCDIR') + '/data/aa-0.csv',
                 glass_file     = os.getenv('PONESRCDIR') + '/data/glass.csv',
                 random_service = None):
        '''
        Define P-OM properties and geometry
        '''
        # ----------------------------------------------------------------------------
        # set up pmt geometry
        # ----------------------------------------------------------------------------
        self.MODULE_RADIUS_M      = 0.2159
        self.PMT_RADIUS_M         = 0.055
        
        # indices in this array correspond to the
        # pmt number in the POM frame -1
        self.PMT_ANGLES = np.array([[57.5, 270.],  # PMT 1
                                    [57.5, 0.],    # PMT 2
                                    [57.5, 90.],   # PMT 3
                                    [57.5, 180.],  # PMT 4
                                    [25., 225.],   # PMT 5
                                    [25., 315.],   # PMT 6
                                    [25., 45.],    # PMT 7
                                    [25., 135.],   # PMT 8
                                    [-57.5, 270.], # PMT 9
                                    [-57.5, 180.], # PMT 10
                                    [-57.5, 90.],  # PMT 11
                                    [-57.5, 0.],   # PMT 12
                                    [-25., 315.],  # PMT 13
                                    [-25., 225.],  # PMT 14
                                    [-25., 135.],  # PMT 15
                                    [-25., 45.]    # PMT 16
                                    ])

        # which PMTs are on top and bottom (hemisphere A, hemisphere B)
        self.TOP_PMTS     = np.ones(len(self.PMT_ANGLES), dtype = bool)
        self.TOP_PMTS[8:] = False

        self.UPPER_RING      = np.ones(len(self.PMT_ANGLES), dtype = bool)
        self.UPPER_RING[4:8] = False
        self.UPPER_RING[12:] = False

        # define home pmts for the k40 characterization
        self.UPPER_HOME_INDEX = 0
        self.LOWER_HOME_INDEX = 4
        self.UPPER_RING_HOME  = self.PMT_ANGLES[self.UPPER_HOME_INDEX]
        self.LOWER_RING_HOME  = self.PMT_ANGLES[self.LOWER_HOME_INDEX]

        # the sum of zenith and azimuth gives a unique pmt identifier.
        # could also maybe use something like two bit base 8 numbers...
        self.ANGLE_SUMS       = np.sum(self.PMT_ANGLES, axis=1)

        # define xyz coordinates for all PMT centres
        x_coordinates = np.multiply(np.sin(np.deg2rad(90. - self.PMT_ANGLES[:, 0])), np.cos(np.deg2rad(self.PMT_ANGLES[:, 1])))
        y_coordinates = np.multiply(np.sin(np.deg2rad(90. - self.PMT_ANGLES[:, 0])), np.sin(np.deg2rad(self.PMT_ANGLES[:, 1])))
        z_coordinates = np.cos(np.deg2rad(90. - self.PMT_ANGLES[:, 0]))

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
    # methods for referencing PMTs
    # ----------------------------------------------------------------------------
    def get_pmt_pair_angle(self, pmt_1, pmt_2):
        '''
        Returns the angle in degrees separating two PMTs
        '''
        angle = np.arccos(self.PMT_MATRIX[pmt_1].dot(self.PMT_MATRIX[pmt_2]))
        angle = np.round(np.rad2deg(angle), 2)

        return angle


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
        aa_data     = np.loadtxt(aa_file, delimiter=',')
        aa_function = interp1d(aa_data.T[0], aa_data.T[1], bounds_error=False, fill_value=0)

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
    

    def get_angular_acceptance(self, pmt_distance_list):
        '''
        Returns an array of angular acceptances
        for a given array of distances between
        pmts and photons
        '''
        distances    = np.hstack(pmt_distance_list)
        efficiencies = np.nan_to_num(self.aa_function(distances))

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
        Returns the PMT that is hit by a given
        photon based on geometry along with the
        distance and angle to that PMT which
        can eventually be used for applying
        angular acceptance
        '''
        photon_position            = photon.pos
        photon_position_normalized = np.array([photon_position.x,photon_position.y, photon_position.z]) / photon_position.magnitude

         # deterine the angle between all the pmt vectors and the photon vector
        pmt_photon_angles = self.PMT_MATRIX.dot(photon_position_normalized)
        
        # check if the ditance at the module radius between the pmt vector
        # and the photon vector falls within the PMT size, if so, mark as a hit
        # return a list of indices that have been 'hit'
        pmt_vector_distances = self.MODULE_RADIUS_M * np.sin(np.arccos(pmt_photon_angles))
        pmts_hit             = np.where(np.logical_and(pmt_vector_distances < self.PMT_RADIUS_M, pmt_photon_angles >= 0))[0]

        if len(pmts_hit) == 1:
            # one pmt hit
            hit_angle    = pmt_photon_angles[pmts_hit]
            hit_distance = pmt_vector_distances[pmts_hit]

            return [pmts_hit[0], hit_distance, hit_angle]
        
        elif len(pmts_hit) > 1:
            # multiple pmts hit
            raise IndexError('Error! Two pmts are being hit at the same time!', pmts_hit)
        
        else:
            # no pmt hit return a PMT label 100 to mark an error
            return [100., 100., 100.]


    def get_probabilities(self, photons, pmt_photon_distances):
        '''
        Returns a list of detection probabilities
        correspon to the given list of photons and
        distances between photons and PMTs
        '''
        weights = np.array([p.weight for p in photons])

        quantum_efficiency    = self.get_quantum_efficiency(photons)
        angular_acceptance    = self.get_angular_acceptance(pmt_photon_distances)
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

        for i, photon in enumerate(photon_list):
            # -------------------------------------------------------------
            # the module definition is rotated 90 degrees (Ti can in the x-y plane)
            # because it's a nicer symmetry to work with this way. however, this
            # means we have to rotate each photon 90 degrees around the x axis
            # to transform it to the coordinate space of the module definition

            # ------ !!!!!! WARNING !!!!!! ------
            # DIRECTION IS NOT CURRENTLY USED IN ACCEPTANCE SO IT IS NOT ROTATED
            # BUT AT SOME POINT IT WILL PROBABLY BE INCLUDED. CAN'T SET DIRECTION XYZ
            # LIKE POSITION BUT NEED TO SET DIRECTION THETA AND PHI INSTEAD

            rotation_matrix   = np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]])
            position          = np.array([photon.pos.x, photon.pos.y, photon.pos.z])
            rotated_position  = np.dot(rotation_matrix, position)

            photon.pos.x = rotated_position[0] * I3Units.m #Should this be x,y,z????? bug?
            photon.pos.x = rotated_position[1] * I3Units.m
            photon.pos.x = rotated_position[2] * I3Units.m
            # -------------------------------------------------------------

            pmt, hit_distance, hit_angle = self.get_pmt(photon)
            
            pmt_list[i]          = pmt
            hit_distance_list[i] = hit_distance
    
        probability_list = self.get_probabilities(photon_list, hit_distance_list)

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
    

    def rotate_pmts_home(self, hit_pmts):
        '''
        Returns an array of PMTs that have been
        rotated/vertically flipped to include the
        home PMT. Used to merge symmetrically
        identical PMT combinations
        '''
        # determine the home pmt we want to use. use the lower
        # ring home if only lower ring pmts were hit. in any other
        # case (some number of upper ring pmts were hit), use the
        # upper ring home
        home                   = self.LOWER_RING_HOME
        upper_ring_pmt_indices = np.where(self.UPPER_RING[hit_pmts])[0]
        if len(upper_ring_pmt_indices != 0):
            # if we need to use an upper ring pmt as the home pmt,
            # move it to position 0 in the array because we will
            # rotate whatever pmt is at position 0 to the home
            # locaion
            home              = self.UPPER_RING_HOME
            first_upper_index = upper_ring_pmt_indices[0]

            hit_pmts[0], hit_pmts[first_upper_index] = hit_pmts[first_upper_index], hit_pmts[0]


        hit_pmt_angles = self.PMT_ANGLES[hit_pmts]
        new_angles     = np.zeros_like(hit_pmt_angles)


        new_angles[:, 0] = hit_pmt_angles[:, 0]
        if not self.TOP_PMTS[hit_pmts[0]]:
            new_angles[:, 0] *= -1
        
        azimuth_difference = hit_pmt_angles[0, 1] - home[1]
        new_angles[:, 1]    = hit_pmt_angles[:, 1] - azimuth_difference

        # make sure the angles are within [0, 360) degrees
        # so we can match with the defined pmt angles
        new_angles[:, 1][np.where(new_angles[:, 1] < 0)[0]] += 360
        new_angles[:, 1][np.where(new_angles[:, 1] >= 360)[0]] -= 360

        transformed_pmts = np.zeros_like(hit_pmts, dtype=int)
        new_angle_sums   = np.sum(new_angles, axis=1)

        for i, angle_sum in enumerate(new_angle_sums):
            transformed_pmts[i] = np.where(self.ANGLE_SUMS == angle_sum)[0][0]
        
        return transformed_pmts
    

    def horizontal_symmetric_pmts(self, hit_pmts):
        '''
        Returns an array of PMTs flipped horizontally
        across the home PMT. Used to further compress
        symmetrically identical PMT combinations
        '''
        hit_pmt_angles = self.PMT_ANGLES[hit_pmts]

        if self.UPPER_HOME_INDEX in hit_pmts:
            home = self.UPPER_RING_HOME
        elif self.LOWER_HOME_INDEX in hit_pmts:
            home = self.LOWER_RING_HOME
        else:
            raise Exception(f'No home PMT found in PMT list! (hit_pmts: {hit_pmts})')

        azimuth_differences = hit_pmt_angles[:, 1] - home[1]
        flipped_azimuths    = hit_pmt_angles[:, 1] - (2 * azimuth_differences)

        hit_pmt_angles[:, 1] = flipped_azimuths

        # make sure the angles are within [0, 360) degrees
        # so we can match with the defined pmt angles
        hit_pmt_angles[:, 1][np.where(hit_pmt_angles[:, 1] < 0)[0]] += 360
        hit_pmt_angles[:, 1][np.where(hit_pmt_angles[:, 1] >= 360)[0]] -= 360

        flipped_pmts       = np.zeros(len(hit_pmts), dtype=int)
        flipped_angle_sums = np.sum(hit_pmt_angles, axis=1)

        for i, angle_sum in enumerate(flipped_angle_sums):
            flipped_pmts[i] = np.where(self.ANGLE_SUMS == angle_sum)[0][0]
        
        return flipped_pmts
