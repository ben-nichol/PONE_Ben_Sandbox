#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from icecube.icetray import I3Units
import FunctionClasses

nominalModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/IceCube/AngularAcceptance.dat', ndmin = 1)
compModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/AngularAcceptance.dat', ndmin = 1)
wavelengthAccNom = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/IceCube/DOMEfficiency.dat', unpack = True)
wavelengthAccComp = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/DOMEfficiency.dat', unpack = True)

nominalModel = FunctionClasses.Polynomial(nominalModel, -1, 1)
compModel = FunctionClasses.Polynomial(compModel, -1, 1)

costheta = np.linspace(-1, 1, 100)

yNom = []
yComp = []

for element in costheta:
    yNom.append(nominalModel.getValue(element))
    yComp.append(compModel.getValue(element))

plt.figure()
plt.plot(costheta, yNom, 'r--', label = "IceCube Model")
plt.plot(costheta, yComp, 'b', label = "MDOM Model")
plt.title("DOM Angular Acceptance")
plt.xlabel(r'$\cos{\beta}$')
plt.ylabel("Acceptance")
plt.legend()

plt.figure()
plt.plot(wavelengthAccNom[0]/I3Units.nanometer, wavelengthAccNom[1], 'r--', label = "IceCube Model")
plt.plot(wavelengthAccComp[0]/I3Units.nanometer, wavelengthAccComp[1], 'b', label = "MDOM Model")
plt.title("DOM Wavelength Acceptance")
plt.xlabel('Wavelength (nm)')
plt.ylabel("Acceptance")
plt.legend()
plt.show()