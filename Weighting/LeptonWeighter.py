from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import LeptonWeighter as LW
from icecube.LeptonInjector import *
import numpy as np
from math import sqrt
from copy import deepcopy
import nuSQuIDS as nsq
import os

class LeptonWeighter(icetray.I3ConditionalModule):
    """
    Wrapper to the opensource lepton weighter.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("crosssectdir","Directory for cross-sections","")
        self.AddParameter("config","Path to the config.lic file","")
        self.AddParameter("fluxConst","Constant for flux parameterization",10**-18)
        self.AddParameter("fluxIndex","Index for flux parameteriazation",1.)
        self.AddParameter("fluxScale","Scale for flux parameterization",10**5)
        self.AddParameter("injectionRadius"," Radius of injection",500.)
        self.AddParameter("computetaueffective","Compue the tau effective weight for muons",False)
        self.AddParameter("nuSQuIDS_WeightTables","Tables for Earth survival weights",os.getenv('PONESRCDIR')+"/data/nsq_allneu_propagation_weight_gamma2.h5")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.crosssectdir = self.GetParameter("crosssectdir")
        self.config = self.GetParameter("config")
        self.fluxConst = self.GetParameter("fluxConst")
        self.fluxIndex = self.GetParameter("fluxIndex")
        self.fluxScale = self.GetParameter("fluxScale")
        self.computetaueffective = self.GetParameter("computetaueffective")

        self.simulation_generators = LW.MakeGeneratorsFromLICFile(self.config)

        self.xs = LW.CrossSectionFromSpline(
                    self.crosssectdir+"/dsdxdy_nu_CC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nubar_CC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nu_NC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nubar_NC_iso.fits")
        self.flux = LW.PowerLawFlux(self.fluxConst, self.fluxIndex, self.fluxScale)
        self.injectionRadius = self.GetParameter("injectionRadius")
        self.weight_event = LW.Weighter( self.flux, self.xs, self.simulation_generators ) 

        self.nsq_atm=nsq.nuSQUIDSAtm(self.GetParameter("nuSQuIDS_WeightTables"))
        self.units=nsq.Const()

        cosrange =  [cth for ic,cth in enumerate(self.nsq_atm.GetCosthRange())]
        self.cosmax = max(cosrange)
        self.cosmin = min(cosrange)
        energyrange = [E for ie,E in enumerate(self.nsq_atm.GetERange())]
        self.Emax = max(energyrange)
        self.Emin = min(energyrange)
        self.units = nsq.Const()
        self.Eflux = lambda E: 1e18*E**(-2.0)

    def Simulation(self,frame) :
        if frame.Has("LeptonInjectorProperties"):
            try:
                self.injectionRadius = frame["LeptonInjectorProperties"].injectionRadius
            except:
                self.injectionRadius = frame["LeptonInjectorProperties"].cylinderRadius

        self.PushFrame(frame)

    def GetnuSQUiDSIndex(self,primary_type) :
        typeindex = 0
        antiindex = 0
        if primary_type == LW.ParticleType.NuMu or primary_type == LW.ParticleType.NuMuBar :
            typeindex =  1
        elif primary_type == LW.ParticleType.NuTau or primary_type == LW.ParticleType.NuTauBar :
            typeindex =  2
        if primary_type == LW.ParticleType.NuEBar or primary_type == LW.ParticleType.NuMuBar or primary_type == LW.ParticleType.NuTauBar :
            antiindex = 1

        return typeindex,antiindex

    def DAQ(self,frame) :

        if not frame.Has("EventProperties"):
            self.PushFrame(frame)

        eventproperties = frame["EventProperties"]
        mctree = frame["I3MCTree"]
         
        neutrino = mctree[0]

        LWevent = LW.Event()
        LWevent.energy = neutrino.energy
        LWevent.zenith = neutrino.dir.zenith
        LWevent.azimuth = neutrino.dir.azimuth
    
        LWevent.interaction_x = eventproperties.finalStateX
        LWevent.interaction_y = eventproperties.finalStateY
        LWevent.final_state_particle_0 = LW.ParticleType(eventproperties.finalType1)
        LWevent.final_state_particle_1 = LW.ParticleType(eventproperties.finalType2)
        LWevent.primary_type = LW.ParticleType(eventproperties.initialType)
        LWevent.radius = self.injectionRadius
        LWevent.total_column_depth = eventproperties.totalColumnDepth
        LWevent.x = neutrino.pos.x
        LWevent.y = neutrino.pos.y
        LWevent.z = neutrino.pos.z

        weight = self.weight_event(LWevent)

        weightone = self.weight_event.get_oneweight(LWevent)

        weight_taueff = 0.0
        weightone_taueff = 0.0

        if self.computetaueffective and (LWevent.primary_type == LW.ParticleType.NuMu or LWevent.primary_type == LW.ParticleType.NuMuBar) :
            weight_taueff = self.get_effective_tau_weight(LWevent)
            weightone_taueff = self.get_effective_tau_oneweight(LWevent)

        survival_weight = 1.0

        nusq_energy = neutrino.energy*self.units.GeV
        cos_zen = np.cos(neutrino.dir.zenith)

        typeindex,antiindex=self.GetnuSQUiDSIndex(LWevent.primary_type)

        if (cos_zen < self.cosmax and cos_zen > self.cosmin) :
            if (nusq_energy < self.Emax and nusq_energy > self.Emin):
                survival_weight = self.nsq_atm.EvalFlavor(typeindex,cos_zen,nusq_energy,antiindex)/self.Eflux(nusq_energy) 
            elif nusq_energy <= self.Emin :
                survival_weight = 1.0
            else :
                survival_weight = 0.0
        elif (nusq_energy < self.Emax and nusq_energy > self.Emin):
            print("wouldn't happen")
            survival_weight = self.nsq_atm.EvalFlavor(typeindex,0.0,nusq_energy,antiindex)/self.Eflux(nusq_energy)
        else :
            #kill event
            survival_weight = 0.0
    
        flux = self.fluxScale*(neutrino.energy*self.units.GeV)**-self.fluxIndex + self.fluxConst

        if weight==np.nan:
            raise ValueError("Bad Weight!")

        frame["LeptonInjection_weight"] = dataclasses.I3Double(weight)
        frame["LeptonInjection_oneweight"] = dataclasses.I3Double(weightone)
        frame["LeptonInjection_weight_taueff"] = dataclasses.I3Double(weight_taueff)
        frame["LeptonInjection_oneweight_taueff"] = dataclasses.I3Double(weightone_taueff)
        frame["LeptonInjection_SurvivalProb"] = dataclasses.I3Double(survival_weight)
        frame["LeptonInjection_flux"] = dataclasses.I3Double(flux)
                
        self.PushFrame(frame)
