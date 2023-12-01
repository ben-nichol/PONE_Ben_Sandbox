import numpy as np 
from icecube.simclasses import I3CLSimFunctionRefIndexQuanFry
from icecube.icetray import I3Tray,I3Units

c = 0.299792458  # speed of light
lambda_s = 120.0  # scattering length of light for violet light
lambda_a = 35  # absorption length of light for violet light
tau = 50.0  # optical attenuation length

indextable = I3CLSimFunctionRefIndexQuanFry(pressure=215.82225*I3Units.bar, temperature=1.78, salinity=34.82*I3Units.perThousand)


def GetIndex(wl=450.) :
    return indextable.GetValue(450.)
    #from http://research.engr.oregonstate.edu/parrish/index-refraction-seawater-and-freshwater-function-wavelength-and-temperature#:~:text=While%20the%20internet%20might%20tell,and%20pressure%20of%20the%20water.
    a = -1.50156e-6
    b = 1.07085e-7
    c = -4.27594e-5
    d = -1.60476e-4
    e = 1.39807
    T = 2.0
    return  a*T**2 + b*wl**2 + c*T + d*wl + e

n = GetIndex()
theta_c = np.arccos(1.0 /n)  # Cherenkov angle in water in radians
c_n = c/n # phase speed
def GetGroupIndex(wl=450.) :
    return indextable.GetValue(450.) - 450.e9*(indextable.GetDerivative(450.))
    a = -1.50156e-6
    b = 1.07085e-7
    c = -4.27594e-5
    d = -1.60476e-4
    e = 1.39807
    T = 2.0
    n = GetIndex(wl)
    dndwl = b*wl*2.0 + d
    return n/(1.0+(wl/n)*dndwl)
ngroup = GetGroupIndex()
c_ngroup = c/ngroup
