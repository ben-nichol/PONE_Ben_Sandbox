from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import LeptonWeighter as LW
from icecube.LeptonInjector import *
import numpy as np
from math import sqrt
from copy import deepcopy

class LeptonWeighter(icetray.I3ConditionalModule):
    """
    Wrapper to the opensource lepton weighter.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("crosssectdir","Directory for cross-sections","")
        self.AddParameter("config","Path to the config.lic file","")
        self.AddParameter("fluxConst","Constant for flux parameterization",10**-18)
        self.AddParameter("fluxIndex","Index for flux parameteriazation",2.)
        self.AddParameter("fluxScale","Scale for flux parameterization",10**5)
        self.AddParameter("injectionRadius"," Radius of injection",500.)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.crosssectdir = self.GetParameter("crosssectdir")
        self.config = self.GetParameter("config")
        self.fluxConst = self.GetParameter("fluxConst")
        self.fluxIndex = self.GetParameter("fluxIndex")
        self.fluxScale = self.GetParameter("fluxScale")

        self.simulation_generators = LW.MakeGeneratorsFromLICFile(self.config)

        self.xs = LW.CrossSectionFromSpline(
                    self.crosssectdir+"/dsdxdy_nu_CC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nubar_CC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nu_NC_iso.fits",
                    self.crosssectdir+"/dsdxdy_nubar_NC_iso.fits")
        self.flux = LW.PowerLawFlux(self.fluxConst, self.fluxIndex, self.fluxScale)
        self.injectionRadius = self.GetParameter("injectionRadius")
        self.weight_event = LW.Weighter( self.flux, self.xs, self.simulation_generators ) 

    def Simulation(self,frame) :
        #if frame.Has("LeptonInjectorProperties"):
        #    self.injectionRadius = frame["LeptonInjectorProperties"].injectionRadius
        self.PushFrame(frame)

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

        if weight==np.nan:
            raise ValueError("Bad Weight!")

        frame["LeptonInjection_weight"] = dataclasses.I3Double(weight)
        frame["LeptonInjection_oneweight"] = dataclasses.I3Double(weightone)            
                
        self.PushFrame(frame)
