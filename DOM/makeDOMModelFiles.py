#!/usr/bin/env python

from icecube import dataclasses
from icecube.icetray import I3Units

import numpy as np
import FunctionClasses
from os.path import expandvars

'''
All functions below taken/modified from clsim code. This modification was needed 
to make the model compatible with new hit detection script.
'''    
    
def GetDOMQuantumEfficiency():
    """
    A function to return the quantum efficiency as instance of
    I3CLSimFunctionFromTable
    """
    # Data copied from the km3 file hit-ini_optic.f
    # BB5912 from Hamamatsu
    q_eff_0 = 0.01
    q_eff_reverse = [1.988, # at 610 nm
                     2.714, # at 600 nm
                     3.496,
                     4.347,
                     5.166,
                     6.004,
                     6.885,
                     8.105,
                     10.13,
                     13.03,
                     15.29,
                     16.37,
                     17.11,
                     17.86,
                     18.95,
                     20.22,
                     21.26,
                     22.10,
                     22.65,
                     23.07,
                     23.14,
                     23.34,
                     22.95,
                     22.95,
                     22.74,
                     23.48,
                     22.59,
                     20.61,
                     17.68,
                     13.18,
                     7.443,
                     2.526  # at 300 nm
                    ]
                    
    q_eff_reverse = [ (q_eff_0 * i) for i in q_eff_reverse ]
    q_eff_reverse.reverse() # reverse the list (in-place)

    startVal = 300.*I3Units.nanometer
    steps = 10.*I3Units.nanometer
    inputs = [(startVal + i*steps) for i in range(len(q_eff_reverse))]

    return FunctionClasses.FunctionFromTable(inputs, q_eff_reverse)


def GetDOMGlassAbsorptionLength():
    """
    A function to return the absoprtion length of
    the glass sphere of an ANTARES OM
    """

    # Data copied from the km3 file hit-ini_optic.f
    # as measured by Pavel
    # the last 3 bins contain measurements from Saclay (6/3/98)
    al_glass_reverse = [148.37, # at 610 nm
                        142.87, # at 600 nm
                        135.64,
                        134.58,
                        138.27,
                        142.40,
                        147.16,
                        151.80,
                        150.88,
                        145.68,
                        139.70,
                        126.55,
                        118.86,
                        113.90,
                        116.08,
                        109.23,
                        81.63,
                        65.66,
                        77.30,
                        73.02,
                        81.25,
                        128.04,
                        61.84,
                        19.23,
                        27.21, 
                        18.09,  
                        8.41,  
                        3.92,  
                        1.82,  
                        0.84,  
                        0.39, 
                        0.17 # at 300 nm
                       ]
                       
    # Apply units
    al_glass_reverse = [ (i * I3Units.cm) for i in al_glass_reverse]
    al_glass_reverse.reverse() # reverse the list (in-place)

    startVal = 300.*I3Units.nanometer
    steps = 10.*I3Units.nanometer
    inputs = [(startVal + i*steps) for i in range(len(al_glass_reverse))]

    return FunctionClasses.FunctionFromTable(inputs, al_glass_reverse)          

def GetDOMGelAbsorptionLength():
    """
    A function to return the absorption length
    the gel of an ANTARES OM
    Note: The file hit-ini_optic.f has three different 
    datasets for this absorption length!
    However in the file hit.f it always is initialized with the
    same (gel_id=1). Thus this one is implemented here.
    """
    
    # Data copied from the km3 file hit-ini_optic.f
    # GEL WACKER (default)
    al_gel_default_reverse = [100.81, # at 610 nm
                              99.94, # at 600 nm
                              99.89, 
                              96.90, 
                              96.42, 
                              94.36, 
                              89.09, 
                              90.10,
                              86.95, 
                              85.88, 
                              84.49, 
                              81.08, 
                              78.18, 
                              76.48, 
                              74.55, 
                              72.31,
                              68.05, 
                              66.91, 
                              64.48, 
                              62.53, 
                              59.38, 
                              56.64, 
                              53.29, 
                              48.96,
                              45.71, 
                              41.88, 
                              37.14, 
                              30.49, 
                              23.08, 
                              15.60,  
                              8.00,  
                              0.0000001 # at 300 nm
                             ]

    # Apply units
    al_gel_default_reverse = [ (i * I3Units.cm) for i in al_gel_default_reverse]
    al_gel_default_reverse.reverse() # reverse the list (in-place)
    
    startVal = 300.*I3Units.nanometer
    steps = 10.*I3Units.nanometer
    inputs = [(startVal + i*steps) for i in range(len(al_gel_default_reverse))]

    return FunctionClasses.FunctionFromTable(inputs, al_gel_default_reverse) 


def GetDOMAcceptance(domRadius = 0.2159*I3Units.m): # 17 inch diameter
    """
    The main function to return the effective area
    of the Antares OM
    """
    # constants taken from km3 (hit-eff_area_pmt.f and hit-transmit.f)
    glass_width = 1.5*I3Units.cm
    gel_width = 1.*I3Units.cm
    pmt_collection_efficiency = 0.9
    pmt_diameter = 9.3 * 0.0254*I3Units.m # 9.3 inch PMT
    
    # Geometrical area of the om profile
    pmt_area = np.pi * (pmt_diameter/2.)**2 #is im square meters
    om_area = np.pi*domRadius**2.
    
    # Get the tables from above
    q_eff = GetDOMQuantumEfficiency()
    abs_glass = GetDOMGlassAbsorptionLength()
    abs_gel = GetDOMGelAbsorptionLength()
    
    # Each of the tables above has 32 bins in the wavelength range of 300nm - 610nm,
    # where the exact value is at the beginning of each bin, means
    # value of bin  0 belongs to 300nm
    # value of bin  1 belongs to 310nm
    # value of bin 31 belongs to 610nm
    # Now combine them
    om_eff_area = [0.] # use a single entry at 290nm to have the same range as other functions(wlen)
    for wavelength in range(300, 611, 10):
        this_abs_glass = abs_glass.getValue(wavelength*I3Units.nanometer)
        this_abs_gel = abs_gel.getValue(wavelength*I3Units.nanometer)
        
        if (this_abs_glass <= 0.) or (this_abs_gel <= 0.):
            current_om_eff_area = 0.
        else:
            current_om_eff_area = pmt_area * \
                                  pmt_collection_efficiency * \
                                  q_eff.getValue(wavelength*I3Units.nanometer) * \
                                  np.exp( -( glass_width / this_abs_glass ) ) * \
                                  np.exp( -( gel_width / this_abs_gel ) )
        
        print current_om_eff_area, om_area
        print pmt_area, pmt_collection_efficiency, q_eff.getValue(wavelength*I3Units.nanometer), np.exp( -( glass_width / this_abs_glass ) ), np.exp( -( gel_width / this_abs_gel ) )

        om_eff_area.append(current_om_eff_area)

    effValues = np.array(om_eff_area)/om_area
    startVal = 300.*I3Units.nanometer
    steps = 10.*I3Units.nanometer
    inputs = [(startVal + i*steps) for i in range(len(effValues))]

    return FunctionClasses.FunctionFromTable(inputs, effValues)

def GetIceCubeDOMAcceptance(domRadius = 0.16510*I3Units.m, efficiency=1.0, highQE=False, coverageFactor = 1):
    """
    this is taken from photonics/lib/efficiency.h:
    
    Adopted into photonics Juli 1 2007 /Johan. Interpolation of data
    from Kotoyo Hoshina <kotoyo.hoshina@icecube.wisc.edu>: 
    Jan 15 K.Hoshina :
    Found a wrong value at 260nm whicn I didn't offer in July 1st 2007.
    Copied a value from 270nm.

          0.0002027029 -> 0.0000064522


    ROMEO wavelength effective area

    This is the table of photo-electron acceptance of the
    IceCube PMT after through the glass+gel + PMT photocathode,
    as a function of wavelength.
    It corresponds to a 0p.e. threshold.
    The injection angle (off-axis angle) is 0deg.
    
    The acceptances are calculated by:
    
    acceptance = NPEs generated by photo-cathode (0P.E threshold)\
              / Nphotons_inject_to_1m^2
    """
    dom2007a_eff_area = [
    0.0000064522,
    0.0000064522,
    0.0000064522,
    0.0000064522,
    0.0000021980,
    0.0001339040,
    0.0005556810,
    0.0016953000,
    0.0035997000,
    0.0061340900,
    0.0074592700,
    0.0090579800,
    0.0099246700,
    0.0105769000,
    0.0110961000,
    0.0114214000,
    0.0114425000,
    0.0111527000,
    0.0108086000,
    0.0104458000,
    0.0099763100,
    0.0093102500,
    0.0087516600,
    0.0083225800,
    0.0079767200,
    0.0075625100,
    0.0066377000,
    0.0053335800,
    0.0043789400,
    0.0037583500,
    0.0033279800,
    0.0029212500,
    0.0025334900,
    0.0021115400,
    0.0017363300,
    0.0013552700,
    0.0010546600,
    0.0007201020,
    0.0004843820,
    0.0002911110,
    0.0001782310,
    0.0001144300,
    0.0000509155]
    dom2007a_eff_area = np.array(dom2007a_eff_area)*I3Units.meter2 # apply units (this is an effective area)
    domArea = np.pi*domRadius**2.
    dom2007a_efficiency = efficiency*(dom2007a_eff_area/domArea)
    
    if highQE:
        wv, rde = np.loadtxt(expandvars('$I3_BUILD/ice-models/resources/models/wavelength/wv.rde')).T
        dom2007a_efficiency *= np.interp(260 + 10*np.arange(len(dom2007a_efficiency)), wv, rde)

    startVal = 260.*I3Units.nanometer
    steps = 10.*I3Units.nanometer
    inputs = [(startVal + i*steps) for i in range(len(dom2007a_efficiency))]

    return  FunctionClasses.FunctionFromTable(inputs, dom2007a_efficiency*coverageFactor)

filePath = '/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/'
filenameAA = 'AngularAcceptance.dat'
filenameDE = 'DOMEfficiency.dat'

aaFile = open(filePath + filenameAA, 'w')
deFile = open(filePath + filenameDE, 'w')

# Can either from icecube paramatarization files, or by hand 
# Polynomial takes care of normalization
# coefficients = np.loadtxt(expandvars('$I3_BUILD/ice-models/resources/models/angsens/as.nominal))
coefficients = np.ndarray([1])

angularAcceptance = FunctionClasses.Polynomial(coefficients, -1, 1)

# code in place to change all aspects of DOM efficiency (glass properties, gel properties, Q.E, 
# PMT coverage) but right now only scaling IceCube efficiency by assumed portions of PMT coverage
# due to lack of data. PMT areas treated as hemi-spheres
#domAcceptance = GetDOMAcceptance()
upgradePMTRadius = 40.25 * I3Units.mm
upgradeDOMRadius = 7.0/12 *I3Units.ft
upgradeNumPMTs = 24
upgradeCoverage = (upgradeNumPMTs * np.pi * upgradePMTRadius**2) / (4 * np.pi * upgradeDOMRadius**2)

currentPMTRadius = 5.0/12 *I3Units.ft
currentDOMRadius = 0.16510*I3Units.m
currentCoverage = (2 * np.pi * currentPMTRadius**2) / (4 * np.pi * currentDOMRadius**2)

print upgradeCoverage / currentCoverage

domAcceptance = GetIceCubeDOMAcceptance(coverageFactor = upgradeCoverage / currentCoverage)

aaFile.write("# This file contains polynomial coefficients for DOM angular Acceptance\n")
for i in range(len(angularAcceptance.coefficients)):
    line = str(angularAcceptance.coefficients[i]) + '\n'
    aaFile.write(line)

aaFile.close()

deFile.write("# This file contains input output values for the DOM efficiency distribution (eff vs. wavelength)\n")
for i in range(len(domAcceptance.inputs)):
    line = str(domAcceptance.inputs[i]) + "\t" + str(domAcceptance.outputs[i]) + "\n"
    deFile.write(line)

deFile.close()
