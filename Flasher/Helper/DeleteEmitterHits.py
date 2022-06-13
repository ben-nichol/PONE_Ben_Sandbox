# Icetray module to initiate an isotropic light source

from icecube import icetray
from icecube.dataclasses import ModuleKey
from icecube.simclasses import I3CompressedPhotonSeriesMap


class DeleteEmitterHits(icetray.I3Module):
    """
    Generate an isotropic photon flash into a DAQ Frame.
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("FlasherKey",
                          "Module key of the flasher position",
                          ModuleKey(1,1))
        self.AddParameter("PhotonSeriesName",
                          "Name of the I3Frame Key the photon flash should be written to",
                          "PhotonSeriesName")
        self.AddOutBox("OutBox")
    
    
    # configuration of the icetray module
    def Configure(self):
        self.flasher_key = self.GetParameter("FlasherKey")
        self.photon_series = self.GetParameter("PhotonSeriesName")
        

    # function to write this into the DAQ frame
    def DAQ(self, frame):
        flasher_pulse_series = frame[self.flasher_pulse_series]
        photon_series = I3CompressedPhotonSeriesMap(frame[self.photon_series])
        
        # pop keys of flasher OMs to remove from output hits
        del photon_series[self.flasher_key]
            
        # delete original hit series
        frame.Delete(self.photon_series)
        
        # insert new photon series
        frame[self.photon_series] = photon_series

        # and push frame
        self.PushFrame(frame, "OutBox")
