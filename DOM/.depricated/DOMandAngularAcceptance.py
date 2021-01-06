import sys
sys.path.insert(0,'/home/users/dhilu/P_ONE_dvirhilu/src')

from icecube import dataclasses, dataio, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from experimentModelCode import FunctionClasses
import matplotlib.pyplot as plt
import numpy as np
from os.path import expandvars
from icecube.clsim import I3CLSimFunctionFromTable


inFolder = '/home/users/akatil/P-ONE/git/PONE_NuTau/DOM/IceCube/'
filenameDomEff = 'DOMEfficiency.dat'
filenameAngAcc = 'AngularAcceptance.dat'

DOMRadius = 0.16510*I3Units.m
DOMOversizeFactor = 5.0
icemodel_efficiency_factor = 1.0
UnshadowedFraction = 1.0
domAcceptance = clsim.GetIceCubeDOMAcceptance(
	                                domRadius = DOMRadius*DOMOversizeFactor,
	                                efficiency=icemodel_efficiency_factor*UnshadowedFraction)

HoleIceParameterization = expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm")
domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=HoleIceParameterization)

angAcc = FunctionClasses.Polynomial(np.loadtxt(inFolder + filenameAngAcc, ndmin = 1), -1, 1)
domEffData = np.loadtxt(inFolder + filenameDomEff, unpack = True)
domEff = FunctionClasses.FunctionFromTable(domEffData[0], domEffData[1])

angle = np.linspace(-1, 1, 10000)
AngAcc = []
for i in angle:
    #AngAcc.append(angAcc.getValue(i))
    AngAcc.append(domAngularSensitivity.GetValue(i))

wavelength = np.linspace(3, 6, 10000)*1e-7
DOMAcc = []
for j in wavelength:
    DOMAcc.append(domAcceptance.GetValue(j))

plt.plot(angle, AngAcc, '.')
plt.xlabel('Angle')
plt.ylabel('Angular Acceptance')
plt.title('Anglular Acceptance of IceCube DOM')
plt.savefig('AngularAccIceCube.pdf', dpi=200)
plt.clf()

plt.plot(wavelength, DOMAcc, '.')
plt.xlabel('Wavelength')
plt.ylabel('Wavelength Acceptance')
plt.title('Wavelength Acceptance of IceCube DOM')
plt.savefig('DOMAccIceCube.pdf', dpi=200)
plt.clf()
