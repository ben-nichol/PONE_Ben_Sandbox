#**************************************
# LeptonInjector modified example script
# orginal by Benjamin Smithers
#**************************************


import LeptonInjector as LI
from math import pi
import os 


#Load cross-sections for pone_offline
xs_folder = os.getenv('PONESRCDIR')+"/CrossSectionModels/csms_differential_v1.0"

# Now, we'll make a new injector for muon tracks 
n_events    = 4
diff_xs     = xs_folder + "/dsdxdy_nu_CC_iso.fits"
total_xs    = xs_folder + "/sigma_nu_CC_iso.fits"
is_ranged   = False #default True
final_1     = LI.Particle.ParticleType.MuMinus
final_2     = LI.Particle.ParticleType.Hadrons
the_injector = LI.Injector( n_events , final_1, final_2, diff_xs, total_xs, is_ranged)



deg = pi/180.

# define some defaults 
minE        = 1000.     # [GeV]
maxE        = 10000.   # [GeV]
gamma       = 2. 
minZenith   = 80.*deg
maxZenith   = 180.*deg
minAzimuth  = 0.*deg
maxAzimuth  = 180.*deg

# construct the controller 
controller  = LI.Controller( the_injector, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)  

# specify the output, earth model
path_to = "/usr/local/LeptonInjector/source/resources/earthparams/"
controller.SetEarthModel("Planet", path_to)
controller.NameOutfile("dataio/data_output.h5")
controller.NameLicFile("dataio/config.lic")

# run the simulation
controller.Execute()
