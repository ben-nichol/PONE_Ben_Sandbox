# ag, 16.06.22.

import numpy as np
from iminuit import Minuit
from icecube import icetray, dataclasses
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.OpticalParameters import c, GetIndex, GetGroupIndex
class TimingFit(icetray.I3Module):

    def gamma_t_only(module, direct, origin):
        s = 1e9; GeV = 1.0; m = 1.0
        n = GetIndex(450.)

        pos = lambda t: dataclasses.I3Position(origin.x + c*t*direct.x,origin.y + c*t*direct.y,origin.z + c*t*direct.z)

        z_close = origin.z
        if (1.0 - direct.z**2) > 0.0 :
            z_close = (   origin.z 
                    - direct.z * (origin.x*direct.x+origin.y*direct.y+origin.z*direct.z) 
                    + direct.z * (module.x * direct.x + module.y * direct.y) 
                  ) / (1.0 - direct.z**2)

        t_close = (  module.x*direct.x + module.y*direct.y + z_close*direct.z - (origin.x*direct.x+origin.y*direct.y+origin.z*direct.z)) / c

        pos_close = pos(t_close)
        d_close = np.sqrt((pos_close.x - module.x)**2 + (pos_close.y - module.y)**2)

        d_gamma = n / np.sqrt(n**2-1) * np.sqrt(d_close**2 + (module.z - z_close)**2 * (1.0-direct.z**2))
        t_gamma = t_close + ((module.z - z_close) * direct.z + (n**2-1)/n * d_gamma) / c
        cos_gamma = (1.0-direct.z**2) * (module.z-z_close)/d_gamma + direct.z/n

        return(t_gamma)

    def module_data(self,photons):
        pulse_time = []
        pulse_pos  = []
        for dom in self.domkeys :
                all_pulses = []
                for pmt in range(self.npmts):
                    omkey = AddPMTKey(dom,pmt+1)
                    if omkey in photons:
                        all_pulses += [p.time for p in photons[omkey]]
                if len(all_pulses) > 0:
                    pulse_time.append(min(all_pulses))
                    pulse_pos.append(self.geomap[omkey].position)
        return pulse_time, pulse_pos

    def quality_function1(t1,t2):
        return (t1-t2)**2
    def quality_function2(t1, t2):
        return 2 * np.sqrt(1+(t1-t2)**2/2) - 2

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('PulseSeriesName',
            'Name of the I3RecoPulseSeriesMap containing the PMT hits',
            "photon_hits_pmt_nonoise")
        self.AddParameter('GeoMapName',
            'Name of the I3OMGeoMap containing the module positions',
            "I3OMGeoMap")
        self.AddParameter('LinefitName',
            'Name of the Linefit result',
            "linefit")
        self.AddParameter('FixVertexRadius','Fix the vertex of the fitter',True)
        self.AddParameter("vertRadius","radius to restrict the vertext to",200.)
        self.AddParameter("output","output name","antares")

    def Configure(self):
        self.pulseseries = self.GetParameter("PulseSeriesName")
        self.geomapname = self.GetParameter("GeoMapName")
        self.seedtrack  = self.GetParameter("LinefitName")
        self.vertexRad = self.GetParameter("vertRadius")
        self.output = self.GetParameter("output")
        self.npass = 0

    def Fit(self,pulse_time,pulse_pos,quality_function,guesstrack) :

        def data_loop(phi, cos_theta, cos_o_theta,o_phi,eventtime):
                qsum = 0
                for t1,pos in zip(pulse_time, pulse_pos):
                    t2 = TimingFit.gamma_t_only(
                        module=pos,
                        direct = dataclasses.I3Position(np.sqrt(1.-cos_theta**2.0)*np.cos(phi),np.sqrt(1.-cos_theta**2.0)*np.sin(phi),cos_theta),
                        origin = dataclasses.I3Position(self.vertexRad*np.sqrt(1.-cos_o_theta**2.0)*np.cos(o_phi),self.vertexRad*np.sqrt(1.-cos_o_theta**2.0)*np.sin(o_phi),self.vertexRad*cos_o_theta)
                        )
                    qsum += quality_function(t1-eventtime,t2)
                return qsum

        data_loop.errordef = Minuit.LEAST_SQUARES
        data_loop.verbose = 1

        m = Minuit(data_loop,
                    phi         = guesstrack.dir.phi,
                    cos_theta       = np.cos(guesstrack.dir.theta),
                    cos_o_theta     = np.cos(guesstrack.pos.theta),
                    o_phi       = guesstrack.pos.phi,
                    eventtime = guesstrack.time)

        m.limits["phi"] = (0, 2*np.pi)
        m.limits["cos_theta"] = (-1.0, 1.0)
        m.limits["o_phi"] = (0, 2*np.pi)
        m.limits["cos_o_theta"] = (-1.0, 1.0)
        m.limits["eventtime"] = (-1000.,2000.)
        m.errors["phi"] = 1.0
        m.errors["cos_theta"] = 0.1
        m.errors["cos_o_theta"] = 0.1
        m.errors["o_phi"] = 1.0
        m.errors["eventtime"] = 50./(0.299/1.35)

        m.migrad()
        res = m.values

        fit       = dataclasses.I3Particle()
        fit.shape = dataclasses.I3Particle.InfiniteTrack
        fit.dir   = dataclasses.I3Direction(
                np.sqrt(1.-res["cos_theta"]**2.0) * np.cos(res["phi"]),
                np.sqrt(1.-res["cos_theta"]**2.0) * np.sin(res["phi"]),
                res["cos_theta"])

        fit.pos   = dataclasses.I3Position(self.vertexRad*np.sqrt(1.-res["cos_o_theta"]**2.0)*np.cos(res["o_phi"]),self.vertexRad*np.sqrt(1.0-res["cos_o_theta"]**2.0)*np.sin(res["o_phi"]),self.vertexRad*res["cos_o_theta"])
        fit.time = res["eventtime"]
        self.npass += 1
        return fit, m.fmin.fval
            

    def Geometry(self,frame) :
        self.geomap  = frame[self.geomapname]
        self.domkeys = set()
        self.npmts = 0
        for dom in self.geomap.keys() :
            self.npmts = max(dom.pmt,self.npmts)
            self.domkeys.add(NoPMTKey(dom))

        maxradius = 0.0
        for dom in self.geomap.keys() :
            pos = self.geomap[dom].position
            radius = np.sqrt(pos.x**2.0+pos.y**2.0+pos.z**2.0)
            maxradius = max(maxradius,radius)

        self.vertexRad += maxradius

        self.PushFrame(frame)

    def Physics(self, frame):
        photons = frame[self.pulseseries]
        DOMKeys = set()
        self.npass = 0
        for dom in photons.keys() :
            domkey = NoPMTKey(dom)
            DOMKeys.add(domkey)

        if len(DOMKeys) < 6:
            return

        guesstrack = None

        for guess, quality_function, name in [
            (self.seedtrack,  TimingFit.quality_function1, self.output+"1"),
            (self.output+"1", TimingFit.quality_function2, self.output+"2")]:
            try:              guesstrack = frame[guess][0]
            except TypeError: guesstrack = frame[guess]

            pulse_time, pulse_pos = self.module_data(photons)

            fit, value = self.Fit(pulse_time,pulse_pos,quality_function,guesstrack)
            frame[name] = fit

        self.PushFrame(frame)
