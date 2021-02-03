from icecube import dataclasses, dataio, icetray, NuFlux
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

def astroFlux(Energy):
    '''
    Values taken from: A measurement of the diffuse astrophysical
    muon neutrino flux using eight years of IceCube data.
    Link: https://inspirehep.net/files/99e2e4c9620a0bb4ddd15ff7749090b4
    '''
    #return (1.01*1e-18)*(Energy/(100*1000))**(-2.19)
    #Seperating ani-neutrino and neutrino flux
    return (0.505*1e-18)*(Energy/(100*1000))**(-2.19)

def simpleWeight(oneWeight, neutrino, lepton):
    '''
    simple Weight is different for CC and NC events
    '''
    flux = NuFlux.makeFlux('honda2006')
    astro_flux = astroFlux(neutrino.energy)

    #neutral currents with NuMu atmospheric flux
    if lepton.type == 16 or lepton.type == 12:
        atm_flux = flux.getFlux(dataclasses.I3Particle.ParticleType.NuMu, neutrino.energy, np.cos(neutrino.dir.zenith))
        return 2*oneWeight*(3*astro_flux + atm_flux)

    #neutral currents with NuMuBar atmospheric flux
    elif lepton.type == -16 or lepton.type == -12:
        atm_flux = flux.getFlux(dataclasses.I3Particle.ParticleType.NuMuBar, neutrino.energy, np.cos(neutrino.dir.zenith))
        return 2*oneWeight*(3*astro_flux + atm_flux)

    #Charged Currents, no atmospheric flux contribution here
    else:
        return 4*oneWeight*astro_flux

def weight(frame, nFiles):
    mctree = frame["I3MCTree"]
    neutrino = frame["NuGPrimary"]
    weightDict = frame["I3MCWeightDict"]
    oneWeight = weightDict['OneWeight']
    primary = mctree.primaries
    lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)
    NEvents = weightDict['NEvents']

    liveTime = 365.24*24*60*60
    sim_weight = simpleWeight(oneWeight, neutrino, lepton)
    eventWeight = sim_weight*liveTime/(NEvents*nFiles)

    return eventWeight

def NEvents():
    '''
    Note: File should be from step 1
    Additional step because I mixed NuE and NuTaus during simulation
    '''
    tau_NEvents = []
    e_NEvents = []
    for i in range(0, 2000):
        readFile = dataio.I3File('/data/p-one/akatil/step_1_medium_water/step_1_'+str(i)+'_PONE_Phase1_NuTau_NuE.i3.gz')
        for frame in readFile:
            mctree = frame["I3MCTree"]
            neutrino = frame["NuGPrimary"]

            if neutrino.type == 16 or neutrino.type == -16:
                tau_NEvents.append(neutrino.type)
            if neutrino.type == 12 or neutrino.type == -12:
                e_NEvents.append(neutrino.type)

    return len(tau_NEvents), len(e_NEvents)
