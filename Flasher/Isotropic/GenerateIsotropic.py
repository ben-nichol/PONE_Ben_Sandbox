# Icetray module to initiate an isotropic light source

import random
import numpy as np

from icecube import icetray
from icecube.icetray import I3Units
from icecube.dataclasses import I3Position, I3Direction, ModuleKey
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

class GenerateIsotropic(icetray.I3Module):
    """
    Generate an isotropic photon flash into a DAQ Frame.
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3Frame Key the photon flash should be written to.",
                          "PhotonFlasherPulseSeries")
        
        self.AddParameter("FlasherKey",
                          "The ModuleKey of the Flasher.",
                          ModuleKey(1,1))
        
        self.AddParameter("FlasherPosition",
                          "The position of the photon source.",
                          I3Position(0,0,0))
        
        self.AddParameter("FlasherPulseType",
                          "The I3CLSimFlasherPulse.FlasherPulseType of the photon flashes.\
                          For a list, see: $PONE_SRC/clsim/public/clsim/I3CLSimFlasherPulse.h",
                          I3CLSimFlasherPulse.FlasherPulseType.Uniform405nm)
        
        self.AddParameter("NumberOfPhotons",
                          "The number of photons to inject from the given position.",
                          1)
        
        self.AddParameter("PulseWidth", 
                          "The pulse FWHM width of the pulse in nanoseconds.",
                          5*I3Units.ns)
        
        self.AddParameter("Seed",
                          "Seed for the random number generator.",
                          1234)
        
        self.AddOutBox("OutBox")
    
    # configuration of the icetray module
    def Configure(self):
        self.series_frame_key = self.GetParameter("FlasherPulseSeriesName")
        self.flasher_key = self.GetParameter("FlasherKey")
        self.flasher_position = self.GetParameter("FlasherPosition")
        self.pulse_type = self.GetParameter("FlasherPulseType")
        self.num_of_photons = self.GetParameter("NumberOfPhotons")
        self.pulse_width = self.GetParameter("PulseWidth")
        self.seed = self.GetParameter("Seed")
        

    # definition of the photon pulse
    def generate_pulse(self, flasher_position, photon_direction, pulse_type,
                       number_of_photons, pulse_width):
        # setup pulse instance
        pulse = I3CLSimFlasherPulse()
        
        # flasher position
        pulse.SetPos(flasher_position)
        
        # initial, arbitrary for isotropic
        pulse.SetDir(photon_direction)
        
        # pulse start time
        pulse.SetTime(0.0*I3Units.ns)
        
        # pulse type
        pulse.SetType(pulse_type)
        
        # number of photons
        pulse.SetNumberOfPhotonsNoBias(number_of_photons)
        
        # pulse FWHM
        pulse.SetPulseWidth(pulse_width)

    	# isotropic in 4pi (full sphere)
        pulse.SetAngularEmissionSigmaPolar(+1.0)
        pulse.SetAngularEmissionSigmaAzimuthal(360. * I3Units.deg)
        
        return pulse


    # function to write this into the DAQ frame
    def DAQ(self, frame):
        random.seed(self.seed)
        pulse_series = I3CLSimFlasherPulseSeries()
        
        # isotropic in 4pi (full sphere)
        photon_direction = I3Direction()
        photon_direction.set_theta_phi(0.,0.)
        pulse = self.generate_pulse(self.flasher_position, photon_direction,
                                    self.pulse_type, self.num_of_photons,
                                    self.pulse_width)
        pulse_series.append(pulse)          
                                
        # push to frame
        frame[self.series_frame_key] = pulse_series
        self.PushFrame(frame, "OutBox")


class GenerateIsotropicHemisphere(icetray.I3Module):
    """
    Generate an isotropic photon flash into a DAQ Frame.
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3Frame Key the photon flash should be written to.",
                          "PhotonFlasherPulseSeries")
        
        self.AddParameter("FlasherKey",
                          "The ModuleKey of the Flasher.",
                          ModuleKey(1,1))
        
        self.AddParameter("FlasherPosition",
                          "The position of the photon source.",
                          I3Position(0,0,0))
        
        self.AddParameter("FlasherPulseType",
                          "The I3CLSimFlasherPulse.FlasherPulseType of the photon flashes.\
                          For a list, see: $PONE_SRC/clsim/public/clsim/I3CLSimFlasherPulse.h",
                          I3CLSimFlasherPulse.FlasherPulseType.Uniform405nm)
        
        self.AddParameter("NumberOfPhotons",
                          "The number of photons to inject from the given position.",
                          1)
        
        self.AddParameter("PulseWidth", 
                          "The pulse FWHM width of the pulse in nanoseconds.",
                          5*I3Units.ns)
        
        self.AddParameter("Seed",
                          "Seed for the random number generator.",
                          1234)
        
        self.AddParameter("Upward",
                          "Flash an upward pointing hemisphere.",
                          1)
        
        self.AddParameter("Downward",
                          "Flash a downward pointing hemisphere.",
                          0)
        
        self.AddOutBox("OutBox")
    
    # configuration of the icetray module
    def Configure(self):
        self.series_frame_key = self.GetParameter("FlasherPulseSeriesName")
        self.flasher_key = self.GetParameter("FlasherKey")
        self.flasher_position = self.GetParameter("FlasherPosition")
        self.pulse_type = self.GetParameter("FlasherPulseType")
        self.num_of_photons = self.GetParameter("NumberOfPhotons")
        self.pulse_width = self.GetParameter("PulseWidth")
        self.upward = bool(self.GetParameter("Upward"))
        self.downward = bool(self.GetParameter("Downward"))
        self.seed = self.GetParameter("Seed")
        

    # definition of the photon pulse
    def generate_pulse(self, flasher_position, photon_direction, pulse_type,
                       number_of_photons, pulse_width):
        # setup pulse instance
        pulse = I3CLSimFlasherPulse()
        
        # flasher position
        pulse.SetPos(flasher_position)
        
        # initial, arbitrary for isotropic
        pulse.SetDir(photon_direction)
        
        # pulse start time
        pulse.SetTime(0.0*I3Units.ns)
        
        # pulse type
        pulse.SetType(pulse_type)
        
        # number of photons
        pulse.SetNumberOfPhotonsNoBias(number_of_photons)
        
        # pulse FWHM
        pulse.SetPulseWidth(pulse_width)

    	# isotropic in 4pi (full sphere)
        pulse.SetAngularEmissionSigmaPolar(+0.0)
        pulse.SetAngularEmissionSigmaAzimuthal(360. * I3Units.deg)
        
        return pulse


    # function to write this into the DAQ frame
    def DAQ(self, frame):
        random.seed(self.seed)
        pulse_series = I3CLSimFlasherPulseSeries()
        
        # isotropic in 2pi (hemisphere)
        up_down = 0 if self.upward else np.pi if self.downward else np.pi/2
        photon_direction = I3Direction()
        photon_direction.set_theta_phi(up_down, 0.)
        pulse = self.generate_pulse(self.flasher_position, photon_direction,
                                    self.pulse_type, self.num_of_photons,
                                    self.pulse_width)
        pulse_series.append(pulse)          
                                
        # push to frame
        frame[self.series_frame_key] = pulse_series
        self.PushFrame(frame, "OutBox")