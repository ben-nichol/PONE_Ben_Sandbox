# Icetray module to initiate an isotropic light source

import math
import random
import numpy as np
from datetime import datetime

from icecube import icetray, dataclasses
from icecube.dataclasses import I3Position, I3Direction
from icecube.icetray import I3Units
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

class GenerateIsotropic(icetray.I3Module):
    """
    Generate an isotropic photon flash into a DAQ Frame.
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("SeriesFrameKey",
                          "Name of the I3Frame Key the photon flash should be written to",
                          "PhotonFlasherPulseSeries")
        self.AddParameter("PhotonPosition",
                          "The position of the photon source.",
                          I3Position(0,0,0))
        self.AddParameter("FlasherPulseType",
                          "The I3CLSimFlasherPulse.FlasherPulseType of the photon flashes.\
                          For a list, see: https://github.com/claudiok/clsim/blob/master/public/clsim/I3CLSimFlasherPulse.h#L59",
                          I3CLSimFlasherPulse.FlasherPulseType.LED405nm)
        self.AddParameter("NumberOfPulses",
                          "The number of pulses to inject from the given position.",
                          1)
        self.AddParameter("NumberOfPhotons",
                          "The number of photons to inject from the given position.",
                          1)
        self.AddParameter("PulseWidth", 
                          "The pulse width of the pulse in nanoseconds",
                          5*I3Units.ns)
        self.AddParameter("Seed",
                          "Seed for the random number generator",
                          1234)
        self.AddParameter("Isotropy",
                          "Using isotropic or hemispheric photon emission",
			  True)
        self.AddParameter("Upward",
                          "Using upward uniform hemispheric photon emission",
			  False)
        self.AddParameter("Downward",
                          "Using downward uniform hemispheric photon emission",
			  False)
        self.AddOutBox("OutBox")
    
    # configuration of the icetray module
    def Configure(self):
        self.series_frame_key = self.GetParameter("SeriesFrameKey")
        self.photon_position = self.GetParameter("PhotonPosition")
        self.pulse_type = self.GetParameter("FlasherPulseType")
        self.num_of_pulses = self.GetParameter("NumberOfPulses")
        self.num_of_photons = self.GetParameter("NumberOfPhotons")
        self.pulse_width = self.GetParameter("PulseWidth")
        self.seed = self.GetParameter("Seed")
        self.isotropy = self.GetParameter("Isotropy")
        self.upward = self.GetParameter("Upward")
        self.downward = self.GetParameter("Downward")


    # definition of the photon pulse
    def generate_pulse(self, photon_position, photon_direction, pulse_type,
                       number_of_photons, pulse_width, isotropy):
        # setup pulse instance
        pulse = I3CLSimFlasherPulse()
        pulse.SetPos(photon_position)
        pulse.SetDir(photon_direction)
        pulse.SetTime(0.0*I3Units.ns)
        pulse.SetType(pulse_type)

        # number of photons
        pulse.SetNumberOfPhotonsNoBias(number_of_photons)
        
        # pulse duration
        pulse.SetPulseWidth(pulse_width)

    	# isotropic in 4pi (full sphere)
        if isotropy:
            pulse.SetAngularEmissionSigmaPolar( 180. * I3Units.deg )
            pulse.SetAngularEmissionSigmaAzimuthal( 360. * I3Units.deg )
        
        # isotropic in 2pi (hemisphere)
        else:
            pulse.SetAngularEmissionSigmaPolar( 90. * I3Units.deg )
            pulse.SetAngularEmissionSigmaAzimuthal( 360. * I3Units.deg )
        return pulse


    # function to write this into the DAQ frame
    def DAQ(self, frame):
        random.seed(self.seed)
        pulse_series = I3CLSimFlasherPulseSeries()
        
        # generate n pulses
        for i in range(self.num_of_pulses):
            
            # isotropic in 4pi (full sphere)
            if self.isotropy:
                photon_direction = I3Direction()
                photon_direction.set_theta_phi(0., 0.) # arbitrary due to isotropy
                pulse = self.generate_pulse(self.photon_position, photon_direction,
                                            self.pulse_type, self.num_of_photons,
                                            self.pulse_width, self.isotropy)
                pulse_series.append(pulse)
    
            # one hemisphere
            else:
                position = self.photon_position
                position = I3Position(*[position[0],
                                        position[1],
                                        position[2]])
                
                # define direction
                direction = I3Direction()
                up_down = 0 if self.upward else np.pi if self.downward else np.nan
                direction.set_theta_phi(up_down, 0.)
                
                # define photon pulses
                pulse = self.generate_pulse(position, direction,
                                            0.5*self.num_of_photons, self.isotropy)
                
                # add pulse
                pulse_series.append(pulse)                         
                                    
        # and push to frame
        frame[self.series_frame_key] = pulse_series
        self.PushFrame(frame, "OutBox")
