# Icetray module to initiate an isotropic light source

import copy

from icecube import icetray
from icecube.clsim import I3CLSimFlasherPulseSeries


class DeleteEmitterHits(icetray.I3Module):
    """
    Generate an isotropic photon flash into a DAQ Frame.
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3Frame Key the pulses are written into",
                          "FlasherPulseSeriesName")
        self.AddParameter("PhotonSeriesName",
                          "Name of the I3Frame Key the photon flash should be written to",
                          "PhotonSeriesName")
        self.AddOutBox("OutBox")
    
    
    # configuration of the icetray module
    def Configure(self):
        self.flasher_pulse_series = self.GetParameter("FlasherPulseSeriesName")
        self.photon_series = self.GetParameter("PhotonSeriesName")
        

    # function to write this into the DAQ frame
    def DAQ(self, frame):
        photon_series = copy.deepcopy(frame[self.photon_series])
        
        # pop keys of flasher OMs to remove from output hits
        for key in self.flasher_pulse_series:
            print(key)
            photon_series.pop(key)
            
        # delete original hit series
        frame.Delete(self.photon_series)
        
        # insert new photon series
        frame[self.photon_series] = photon_series

        # and push frame
        self.PushFrame(frame, "OutBox")
