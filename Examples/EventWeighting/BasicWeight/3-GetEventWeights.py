import LeptonWeighter as LW
import h5py as h5
import numpy as np
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Calculate event weights with configurable input files.")
parser.add_argument('--lic', default='mcgen.lic', help='Path to LIC config file (default: mcgen.lic)')
parser.add_argument('--cross-section-location', default='/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/', help='Directory containing cross section files')
parser.add_argument('--input', default='mcgen.h5', help='Input HDF5 file (default: mcgen.h5)')
args = parser.parse_args()

# Define flux parameters
#                         GeV      unitless        GeV
flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   = 3.1536e7 #s

# Create generator
net_generation = LW.MakeGeneratorsFromLICFile(args.lic)

# Create flux, load cross sections 
# This cross section object takes four differential cross sections (dS/dEdxdy) 
#            Neutrino CC-DIS xs
#       Anti-Neutrino CC-DIS xs
#            Neutrino NC-DIS xs
#       Anti-Neutrino NC-DIS xs
cross_section_location = args.cross_section_location
xs = LW.CrossSectionFromSpline(
    cross_section_location + "/dsdxdy_nu_CC_iso.fits",
    cross_section_location + "/dsdxdy_nubar_CC_iso.fits",
    cross_section_location + "/dsdxdy_nu_NC_iso.fits",
    cross_section_location + "/dsdxdy_nubar_NC_iso.fits"
)
flux = LW.PowerLawFlux(flux_params['constant'], flux_params['index'], flux_params['scale'])

# build weighter
weight_event = LW.Weighter(flux, xs, net_generation)

def get_weight(props):
    """
    Accepts the properties list of an event and returns the weight
    """
    LWevent = LW.Event()
    LWevent.energy = props[5]
    LWevent.zenith = props[6]
    LWevent.azimuth = props[7]
    
    LWevent.interaction_x = props[8]
    LWevent.interaction_y = props[9]
    LWevent.final_state_particle_0 = LW.ParticleType(props[10])
    LWevent.final_state_particle_1 = LW.ParticleType(props[11])
    LWevent.primary_type = LW.ParticleType(props[12])
    LWevent.radius = 0 #?
    LWevent.total_column_depth = props[16]
    LWevent.x = 0
    LWevent.y = 0
    LWevent.z = 0
    
    weight = weight_event(LWevent)
    if weight == np.nan:
        raise ValueError("Bad Weight!")
    return weight * liveTime

# load data
data_file = h5.File(args.input)
for event in range(len(data_file['EventProperties'])):
    print("Event Weight: {}".format(get_weight(data_file['EventProperties'][event])))
data_file.close()
