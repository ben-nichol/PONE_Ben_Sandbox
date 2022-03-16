import numpy as np 

c = 0.299792458  # speed of light
n = 1.34  # index of refraction
c_n = c/n # phase speed
ngroup = 1.35557  # group index of refraction
c_ngroup = c/ngroup # group velocity
theta_c = np.arccos(1.0 / n)  # Cherenkov angle in water in radians
lambda_s = 120.0  # scattering length of light for violet light
lambda_a = 35  # absorption length of light for violet light
tau = 50.0  # optical attenuation length
