# ag, 13.06.22.

from os import getenv
import numpy as nn
import pandas as pd
from scipy import integrate, interpolate
from scipy.spatial.transform import Rotation as R

from icecube import icetray, dataclasses
from icecube.icetray import I3Units


class MuonGenerator(icetray.I3Module):
    "Creates atmospheric muons with energy/angle distributions from a CSV file. For each frame, one muon is generated."

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter(
            "Fluxfile",
            "CSV file with flux values from MUTE.",
            getenv("PONESRCDIR") + "/data/AtmosphericMuons/fluxtable_1460m.csv",
        )
        self.AddParameter(
            "InjectionHeight",
            "Height at which muons are created. Should be at least 200m above the detector. Must match with the depth for which the fluxfile was generated",
            700 * I3Units.m,
        )
        self.AddParameter(
            "InjectionRadius",
            "Radius of a sphere that encompasses the detector, with at least 200m clearance to all sides.",
            700 * I3Units.m,
        )
        self.AddParameter("OutputPrefix", "", "MuonGenerator")
        self.AddParameter("RandomService", "Random number generator", None)

    def Configure(self):
        # csv file uses second, cm, MeV
        d = pd.read_csv(
            self.GetParameter("Fluxfile"),
            comment="#",
            delim_whitespace=True,
            header=None,
            nrows=None,
        )

        flux_angle_rad = nn.array(d.loc[0, 1:], dtype="float") / 180 * nn.pi
        flux_energy_mev = nn.array(d.loc[1:, 0], dtype="float")
        flux = nn.array(d.loc[1:, 1:], dtype="float")
        flux[flux < 0] = 0

        angle_int = integrate.simps(flux, flux_energy_mev, axis=0)
        energy_int = integrate.simps(flux, flux_angle_rad, axis=1)

        angle_cflux = nn.cumsum(angle_int * nn.sin(flux_angle_rad)[None, :])
        angle_cflux /= nn.max(angle_cflux)

        energy_cflux = nn.cumsum(flux * flux_energy_mev[:, None], axis=0)
        norm = nn.max(energy_cflux, axis=0)
        norm[norm == 0] = 1
        energy_cflux /= norm
        ang, E = nn.meshgrid(flux_angle_rad, flux_energy_mev)

        self.total_flux = (
            integrate.simps(
                angle_int * nn.sin(flux_angle_rad) * 2 * nn.pi, flux_angle_rad, axis=0
            )
            / I3Units.s
            / I3Units.cm2
        )
        self.angle_interpolator = interpolate.interp1d(
            angle_cflux, flux_angle_rad, kind="linear"
        )
        self.energy_interpolator = interpolate.LinearNDInterpolator(
            (ang.reshape((-1)), energy_cflux.reshape((-1))),
            E.reshape((-1)),
            rescale=True,
        )
        self.flux_interpolator = interpolate.LinearNDInterpolator(
            (ang.reshape((-1)), E.reshape((-1))), flux.reshape((-1)), rescale=True
        )
        self.has_sframe = False

    def DAQ(self, frame):
        rand = self.GetParameter("RandomService")
        radius = self.GetParameter("InjectionRadius")
        height = self.GetParameter("InjectionHeight")
        prefix = self.GetParameter("OutputPrefix")

        angle = self.angle_interpolator(rand.uniform(0, 1))
        energy = self.energy_interpolator(angle, rand.uniform(0, 1))
        flux = float(self.flux_interpolator(angle, energy)) / I3Units.s / I3Units.cm2

        direction = nn.array([0, 0, -1])
        origin = nn.array([nn.sqrt(rand.uniform(0, radius**2)), 0, 0])
        for axis, angle in [
            ("z", rand.uniform(0, 2 * nn.pi)),
            ("y", angle),
            ("z", rand.uniform(0, 2 * nn.pi)),
        ]:
            rot = R.from_euler(axis, angle).as_dcm()
            origin = nn.einsum("bc,c->b", rot, origin)
            direction = nn.einsum("bc,c->b", rot, direction)
        scale = (height - origin[2]) / direction[2]
        origin += scale * direction

        particle = dataclasses.I3Particle()
        particle.type = dataclasses.I3Particle.MuMinus
        particle.energy = energy * I3Units.MeV
        particle.pos = dataclasses.I3Position(*origin)
        particle.dir = dataclasses.I3Direction(*direction)
        particle.time = 0 * I3Units.ns
        particle.length = float("nan")
        particle.shape = dataclasses.I3Particle.ParticleShape.MCTrack
        particle.location_type_string = "InIce"

        tree = dataclasses.I3MCTree()
        tree.add_primary(particle)

        if not self.has_sframe:
            sframe = icetray.I3Frame("S")
            sframe.Put(prefix + "TotalFlux", dataclasses.I3Double(self.total_flux))
            sframe.Put(
                prefix + "Area",
                dataclasses.I3Double(nn.pi * self.GetParameter("InjectionRadius") ** 2),
            )
            self.PushFrame(sframe)
            self.has_sframe = True

        frame.Put(prefix + "I3MCTree", tree)
        frame.Put(prefix + "IndividualFlux", dataclasses.I3Double(flux))
        self.PushFrame(frame)


class MuonReweighter(MuonGenerator):
    "Reweights muons created with MuonGenerator. Takes CSV file with new flux, and adds I3Double to each frame with new muon weight."

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter(
            "Fluxfile",
            "CSV file with flux values from MUTE.",
            getenv("PONESRCDIR") + "/data/AtmosphericMuons/fluxtable_1460m.csv",
        )
        self.AddParameter(
            "InputI3MCTreeName",
            "Name of the I3MCTree in which the muons are stored.",
            "MuonGeneratorI3MCTree",
        )
        self.AddParameter(
            "InputFluxName",
            "Flux of the individual muons.",
            "MuonGeneratorIndividualFlux",
        )
        self.AddParameter("OutputPrefix", "", "MuonReweighter")

    def Configure(self):
        MuonGenerator.Configure(self)

    def Simulation(self, frame):
        prefix = self.GetParameter("OutputPrefix")
        frame.Put(prefix + "TotalFlux", dataclasses.I3Double(self.total_flux))
        self.PushFrame(frame)

    def DAQ(self, frame):
        assert (
            len(frame[self.GetParameter("InputI3MCTreeName")]) == 1
        ), "There should be exactly one muon in the I3MCTree"

        muon = frame[self.GetParameter("InputI3MCTreeName")][0]
        angle = nn.arccos(-muon.dir.z)
        energy = muon.energy / I3Units.MeV
        flux = frame[self.GetParameter("InputFluxName")].value / I3Units.s / I3Units.cm2
        weight = float(self.flux_interpolator(angle, energy)) / flux
        prefix = self.GetParameter("OutputPrefix")

        frame.Put(prefix + "Weight", dataclasses.I3Double(weight))
        self.PushFrame(frame)
