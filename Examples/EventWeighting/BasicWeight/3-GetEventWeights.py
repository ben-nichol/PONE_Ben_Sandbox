import LeptonWeighter as LW
import h5py as h5
import numpy as np
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Calculate event weights with multiple LIC and HDF5 input files.")
parser.add_argument('--lic', nargs='+', required=True, help='One or more LIC config files')
parser.add_argument('--cross-section-location', default='/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/', help='Directory containing cross section files')
parser.add_argument('--input', nargs='+', required=True, help='One or more input HDF5 files')
args = parser.parse_args()

# Define flux parameters
#                         GeV      unitless        GeV
flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   = 3.1536e7 #s


xs = LW.CrossSectionFromSpline(
    args.cross_section_location + "/dsdxdy_nu_CC_iso.fits",
    args.cross_section_location + "/dsdxdy_nubar_CC_iso.fits",
    args.cross_section_location + "/dsdxdy_nu_NC_iso.fits",
    args.cross_section_location + "/dsdxdy_nubar_NC_iso.fits"
)
flux = LW.PowerLawFlux(flux_params['constant'], flux_params['index'], flux_params['scale'])

def get_weight(props, weight_event):
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

#get lic details, each event is weighted against all generators
generators = []
for lic_file in args.lic:
    generators.extend(LW.MakeGeneratorsFromLICFile(lic_file))

#get combined weighter
weight_event = LW.Weighter(flux, xs, generators)

#loop over input files and get weights 
for h5_file in args.input:
    print(f"  Processing HDF5 file: {h5_file}")
    data_file = h5.File(h5_file)
    for event in range(len(data_file['EventProperties'])):
        print("Event Weight: {}".format(get_weight(data_file['EventProperties'][event], weight_event)))
    data_file.close()
