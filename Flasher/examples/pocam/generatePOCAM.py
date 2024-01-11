# Icetray module to initiate a POCAM as a light source

from icecube import icetray, dataclasses
from I3Tray import I3Units
from icecube.dataclasses import *
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
from datetime import datetime
import math
import random
import numpy as np


class GeneratePOCAM(icetray.I3Module):
    """
    Generate a POCAM Flash into a DAQ Frame.
    """

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter(
            "SeriesFrameKey",
            "Name of the I3Frame Key the photon flash should be written to",
            "PhotonFlasherPulseSeries",
        )
        self.AddParameter(
            "PhotonPosition", "The position of the photon source.", I3Position(0, 0, 0)
        )
        self.AddParameter(
            "NumberOfPhotons",
            "The number of photons to inject from the given position.",
            1,
        )
        self.AddParameter(
            "PulseWidth", "The pulse width of the pulse in nanoseconds", 5 * I3Units.ns
        )
        self.AddParameter(
            "FlasherPulseType",
            "The I3CLSimFlasherPulse.FlasherPulseType of the photon flashes. For a list, see: https://github.com/claudiok/clsim/blob/master/public/clsim/I3CLSimFlasherPulse.h#L59",
            I3CLSimFlasherPulse.FlasherPulseType.LED405nm,
        )
        self.AddParameter("Seed", "Seed for the random number generator", 1234)
        self.AddParameter(
            "Isotropy", "Using isotropic or hemispheric photon emission", True
        )
        self.AddOutBox("OutBox")

    # configuration of the icetray module
    def Configure(self):
        self.series_frame_key = self.GetParameter("SeriesFrameKey")
        self.photon_position = self.GetParameter("PhotonPosition")
        self.num_of_photons = self.GetParameter("NumberOfPhotons")
        self.pulse_width = self.GetParameter("PulseWidth")
        self.pulse_type = self.GetParameter("FlasherPulseType")
        self.seed = self.GetParameter("Seed")
        self.isotropy = self.GetParameter("Isotropy")

    # definition of the photon pulse
    def generate_pulse(
        self,
        photon_position,
        photon_direction,
        number_of_photons,
        pulse_width,
        isotropy,
    ):
        pulse = I3CLSimFlasherPulse()
        pulse.SetPos(photon_position)
        pulse.SetDir(photon_direction)
        pulse.SetTime(0.0 * I3Units.ns)
        pulse.SetNumberOfPhotonsNoBias(number_of_photons)
        pulse.SetType(self.pulse_type)

        # pulse duration, 5ns for POCAMs
        pulse.SetPulseWidth(pulse_width)

        # isotropic in 4pi (full sphere)
        if isotropy:
            pulse.SetAngularEmissionSigmaPolar(180.0 * I3Units.deg)
            pulse.SetAngularEmissionSigmaAzimuthal(360.0 * I3Units.deg)

        # isotropic in 2pi (hemisphere)
        else:
            pulse.SetAngularEmissionSigmaPolar(90.0 * I3Units.deg)
            pulse.SetAngularEmissionSigmaAzimuthal(360.0 * I3Units.deg)
        return pulse

    # function to write this into the DAQ frame
    def DAQ(self, frame):
        random.seed(self.seed)
        pulse_series = I3CLSimFlasherPulseSeries()

        # isotropic in 4pi (full sphere)
        if self.isotropy:
            photon_direction = I3Direction()
            photon_direction.set_theta_phi(
                0.0, 0.0
            )  # direction arbitrary due to isotropy
            pulse = self.generate_pulse(
                self.photon_position,
                photon_direction,
                self.num_of_photons,
                self.pulse_width,
                self.isotropy,
            )
            pulse_series.append(pulse)

        # two hemispheres, both isotropic in 2pi, separated by some distance
        # this resembles a realistic pocam instrument
        else:
            pocam_position = self.photon_position
            pocam_position1 = I3Position(
                *[pocam_position[0], pocam_position[1], pocam_position[2] + 0.125]
            )
            pocam_position2 = I3Position(
                *[pocam_position[0], pocam_position[1], pocam_position[2] - 0.125]
            )

            # define emission directions of two hemispheres:
            photon_direction1 = I3Direction()
            photon_direction2 = I3Direction()
            photon_direction1.set_theta_phi(0.0, 0.0)  # one hemisphere points upwards
            photon_direction2.set_theta_phi(np.pi, 0.0)  # one downwards

            # define photon pulses for both
            pulse1 = self.generate_pulse(
                pocam_position1,
                photon_direction1,
                0.5 * self.num_of_photons,
                self.isotropy,
            )
            pulse2 = self.generate_pulse(
                pocam_position2,
                photon_direction2,
                0.5 * self.num_of_photons,
                self.isotropy,
            )
            pulse_series.append(pulse1)
            pulse_series.append(pulse2)

        # and push to frame
        frame[self.series_frame_key] = pulse_series
        self.PushFrame(frame, "OutBox")
