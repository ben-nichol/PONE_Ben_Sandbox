import LeptonWeighter as LW
import h5py as h5
import numpy as np

"""
A Basic Weighting example 
This calculates the weight of each event in the standalone LeptonInjector example script. 

Ben Smithers
benjamin.smithers@mavs.uta.edu

To convert this to one working with i3-LeptonInjector, you will need to modify the event loop and the get_weight function. 
"""

#                         GeV      unitless        GeV
flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   = 3.1536e7 #s

# Create generator
#    if there were multiple LIC files, you would instead make a list of Generators
net_generation = LW.MakeGeneratorsFromLICFile("config.lic")

# Create flux, load cross sections 
# This cross section object takes four differential cross sections (dS/dEdxdy) 
#            Neutrino CC-DIS xs
#       Anti-Neutrino CC-DIS xs
#            Neutrino NC-DIS xs
#       Anti-Neutrino NC-DIS xs
cross_section_location = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/"
xs = LW.CrossSectionFromSpline(
                    cross_section_location+"/dsdxdy_nu_CC_iso.fits",
                    cross_section_location+"/dsdxdy_nubar_CC_iso.fits",
                    cross_section_location+"/dsdxdy_nu_NC_iso.fits",
                    cross_section_location+"/dsdxdy_nubar_NC_iso.fits")
flux = LW.PowerLawFlux( flux_params['constant'] , flux_params['index'] , flux_params['scale'] )

# build weighter
weight_event = LW.Weighter( flux, xs, net_generation )

# Write utility function to get the weight 
def get_weight( props ):
    """
    Accepts the properties list of an event and returns the weight
    """
    LWevent = LW.Event()
    LWevent.energy = props[5]
    LWevent.zenith = props[6]
    LWevent.azimuth = props[7]
    
    LWevent.interaction_x = props[8]
    LWevent.interaction_y = props[9]
    LWevent.final_state_particle_0 = LW.ParticleType( props[10] )
    LWevent.final_state_particle_1 = LW.ParticleType( props[11] )
    LWevent.primary_type = LW.ParticleType( props[12] )
    LWevent.radius = 0 #?
    LWevent.total_column_depth = props[16]
    LWevent.x = 0
    LWevent.y = 0
    LWevent.z = 0
    
    weight = weight_event(LWevent)

    # this would alert us that something bad is happening 
    if weight==np.nan:
        raise ValueError("Bad Weight!")

    return( weight*liveTime )


# load data
data_file = h5.File("output.h5")

for event in range(len( data_file['EventProperties'] )):
    print("Event Weight: {}".format( get_weight( data_file['EventProperties'][event]) ))

data_file.close()
