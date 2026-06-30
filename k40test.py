#!/usr/bin/python

'''
Script that generates an I3 file of input particles
'''


import numpy as np
import scipy.stats as st
from icecube import dataclasses, phys_services, dataio, icetray
from icecube.icetray import I3Units
import os, sys
from optparse import OptionParser



usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option('-o', '--output_file',default='./test_output.i3.gz',
                    dest='OUTPUT_FILE',help='Write output to OUTFILE (.i3{.gz} format')
parser.add_option('-f', '--num_frames', type='int', default='2', dest='NUM_FRAMES',
                    help='Number of frames/MCTrees you would like to make')
parser.add_option('-r', '--radius', type='float', default='50.0', dest='RADIUS',
                    help='World radius to generate electrons in')
parser.add_option('-t', '--frame_length', type='float', default='1.0', dest='FRAME_LENGTH_MS',
                    help='Time length of 1 frame in milliseconds, default 1ms')
parser.add_option('-d', '--depth_offset', type='float', default='0.0', dest='DEPTH',
                    help='Static depth offset to add to all the generated particle positions')

(options,args) = parser.parse_args()



def GetMeanEventRate(salinity = 3.482):
    '''
    Returns the expected K40 decay rate as a function of
    ocean salinity %.

    The returned rate is in [decays / ms / m3]
    '''
    #---------------------------------------------
    # Define constants
    rk  = 1.11         # Potassium fraction in sea salt [%]
    ri  = 0.0117       # K40 isotope fraction [%]
    p   = 1.013        # Ocean water density [g/cm3]
    NA  = 6.022e23     # Avogadro's Number [/mol]
    A   = 39.96        # K40 atomic weight [g/mol]
    hl  = 1.251e9      # K40 half life [years]

    #---------------------------------------------
    # Calculate
    numerator   = (salinity/100) * (rk/100) * (ri/100) * p * (1e6) * NA * np.log(2)
    denominator = A * hl * (365 * 24 * 60 * 60 * 1000)
    mean_rate   = (numerator / denominator)

    return mean_rate



def GetElectronEnergies(number_of_electrons):
    '''
    Returns an array of energies [MeV] sampled from
    the decay electron energy spectrum
    '''
    electron_data = np.genfromtxt('/home/jakubs/projects/def-mdanning/jakubs/k40/simulation/input/scripts/resources/electron-energy-distribution.csv', delimiter=',')
    energies      = electron_data[:, 0]
    probabilities = electron_data[:, 1]

    return np.random.choice(energies, size=number_of_electrons, p=probabilities)



def GetDecayTimes(frame_length_ms, event_rate, expected_total_events):
    '''
    Returns an array of times [ns] of decay events
    within the given frame length [ms] according to the
    given event rate and expected total number of events.
    '''
    time_offsets = np.random.gamma(1, 1/event_rate, expected_total_events * 2)
    time_elapsed = np.cumsum(time_offsets)
    cutoff_index = np.where(time_elapsed > frame_length_ms)[0]

    while len(cutoff_index) < 1:
        time_offsets = np.concatenate((time_offsets, np.random.gamma(1, 1/event_rate, expected_total_events * 2)))
        time_elapsed = np.cumsum(time_offsets)
        cutoff_index = np.where(time_elapsed > frame_length_ms)[0]

    return time_elapsed[0:cutoff_index[0]] * 1e6 # convert ms to ns



def GetDecayPositions(number_of_decays, radius=options.RADIUS):
    '''
    Returns an array of xyz coordinates corresponding to
    decay positions. There are sampled uniformly through
    a spherical volume with the given radius [m].
    '''
    random_coordinates   = np.random.normal(0., 1., 3 * number_of_decays).reshape(number_of_decays, 3)
    scaled_radius        = radius * np.random.rand(number_of_decays)**(1.0 / 3.0)
    radius_normalization = np.linalg.norm(random_coordinates, axis=1)

    particle_positions        = scaled_radius[:, np.newaxis] * random_coordinates / radius_normalization[:, np.newaxis]
    particle_positions[:, 2] += options.DEPTH

    return particle_positions



def GetDecayDirection(number_of_decays):
    '''
    Returns an array of xyz coordinates corresponding to
    decay direction unit vectors sampled uniformly in 3 space.
    '''
    random_coordinates   = np.random.normal(0., 1., 3 * number_of_decays).reshape(number_of_decays, 3)
    radius_normalization = np.linalg.norm(random_coordinates, axis=1)

    return random_coordinates / radius_normalization[:, np.newaxis]



# Define constants

# we will only account for the fraction of all k40 decays that go
# through beta decay or electron capture
BR_BETA = 0.8928 # branching ratio for beta decaay
BR_EC   = 0.1067 # branching ratio for electron capture

EC_FRACTION = BR_EC / (BR_BETA + BR_EC)

GAMMA_ENERGY_MEV        = 1.46083
CHERENKOV_THRESHOLD_MEV = 0.7


decay_rate_m_ms = GetMeanEventRate()
world_volume_m  = (4. / 3.) * np.pi * options.RADIUS**3
event_rate_ms   = decay_rate_m_ms * world_volume_m * (BR_BETA + BR_EC)

expected_frame_events = int(decay_rate_m_ms * options.FRAME_LENGTH_MS * world_volume_m)

with dataio.I3File(options.OUTPUT_FILE,'w') as output_file:
    for f in np.arange(options.NUM_FRAMES):
        frame = icetray.I3Frame('Q')
        tree  = dataclasses.I3MCTree()


        # ------------------------------------------------------------------------------------
        # first make a dummy primary particle, all the decay products that
        # need to be propagated will be children secondaries of this particle
        # make a filler primary particle
        primary_particle        = dataclasses.I3Particle()
        primary_particle.type   = dataclasses.I3Particle.EMinus
        primary_particle.energy = 0.1 * I3Units.MeV
        primary_particle.pos    = dataclasses.I3Position(100. * I3Units.m, 100. * I3Units.m, 100. * I3Units.m)
        primary_particle.dir    = dataclasses.I3Direction(1., 1., 1.)
        primary_particle.time   = 0. * I3Units.ns
        primary_particle.length = float('nan')

        primary_particle.location_type_string = 'InIce'
        primary_particle.shape_string         = 'Primary'
        tree.add_primary(primary_particle)
        # ------------------------------------------------------------------------------------


        decay_times = GetDecayTimes(options.FRAME_LENGTH_MS, event_rate_ms, expected_frame_events)
        num_decays  = len(decay_times)

        # filter any beta decays that have energies below the cherenkov threshold
        electron_energies = GetElectronEnergies(num_decays)
        ec_decay_indices  = np.where(np.random.rand(num_decays) < EC_FRACTION)[0]

        beta_mask                   = np.ones(electron_energies.size, dtype=bool)
        beta_mask[ec_decay_indices] = False
        beta_decay_indices          = np.arange(len(electron_energies))[beta_mask]

        sub_cherenkov_indices                     = np.where(electron_energies < CHERENKOV_THRESHOLD_MEV)[0]
        sub_cherenkov_mask                        = np.zeros(electron_energies.size, dtype=bool)
        sub_cherenkov_mask[sub_cherenkov_indices] = True

        sub_cherenkov_electron_mask = np.logical_and(beta_mask, sub_cherenkov_mask)
        filtered_indices            = np.invert(sub_cherenkov_electron_mask)

        decay_times       = decay_times[filtered_indices]
        electron_energies = electron_energies[filtered_indices]
        num_decays        = len(decay_times)
        decay_positions   = GetDecayPositions(num_decays)
        decay_directions  = GetDecayDirection(num_decays)

        for d in np.arange(num_decays):
            if d in ec_decay_indices:
                particle_type   = dataclasses.I3Particle.Gamma
                particle_energy = GAMMA_ENERGY_MEV
            else:
                particle_type   = dataclasses.I3Particle.EMinus
                particle_energy = electron_energies[d]

            particle        = dataclasses.I3Particle()
            particle.type   = particle_type
            particle.energy = particle_energy * I3Units.MeV
            particle.pos    = dataclasses.I3Position(decay_positions[d][0] * I3Units.m, decay_positions[d][1] * I3Units.m, decay_positions[d][2] * I3Units.m)
            particle.dir    = dataclasses.I3Direction(decay_directions[d][0], decay_directions[d][1], decay_directions[d][2])
            particle.time   = decay_times[d] * I3Units.ns
            particle.length = float('nan')

            particle.location_type_string = 'InIce'
            particle.shape_string         = 'Cascade'
            tree.append_child(primary_particle, particle)

        frame.Put('I3DecayProducts', tree)
        output_file.push(frame)