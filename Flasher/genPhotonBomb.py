from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
import numpy as np

from I3Tray import I3Units

class PhotonBomb(icetray.I3ConditionalModule):
    """
    Generates I3CLSimFlasherPulse objects for 370nm and 405nm LEDs (I&II)

    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3CLSimFlasherPulseSeries to write",
                          "I3CLSimFlasherPulseSeries")
        self.AddParameter("PhotonsPerPulse",
                          "Photons to emit per pulse",
                          5)
        self.AddParameter("NumPulses",
                          "Num Pulses",
                          5)
        self.AddParameter("FlashTime",
                          "Time (within each event) at which to flash",
                          0.*I3Units.ns)
        self.AddParameter("RandomService","Random Service")
        self.AddParameter("Position","Position of photon bomb")
        self.AddParameter("Radius","Radius for simulation",200.)
       
        self.AddOutBox("OutBox")

    def Configure(self):
        self.flasherPulseSeriesName = self.GetParameter("FlasherPulseSeriesName")
        self.photonsPerPulse = self.GetParameter("PhotonsPerPulse")
        self.flashTime = self.GetParameter("FlashTime")
        self.randomService = self.GetParameter("RandomService")
        self.position = self.GetParameter("Position")
        self.numPulses = self.GetParameter("NumPulses")
        self.radius = self.GetParameter("Radius")

    def DAQ(self, frame):
        outputSeries = I3CLSimFlasherPulseSeries()

        numPhotons = self.photonsPerPulse

        for i in range(self.numPulses) :
            newPulse = I3CLSimFlasherPulse()
            newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
            theta = random_service.uniform(0.,np.pi)
            phi = random_service.uniform(0.,2.*np.pi)
            newPulse.pos = dataclasses.I3Position(self.radius*np.cos(phi)*np.sin(theta)*I3Units.m, 
                                                  self.radius*np.sin(phi)*np.sin(theta)*I3Units.m, 
                                                  self.radius*np.cos(theta)*I3Units.m)
            theta = random_service.uniform(0.,np.pi)
            phi = random_service.uniform(0.,2.*np.pi)
            newPulse.dir = dataclasses.I3Direction(np.cos(phi)*np.sin(theta)*I3Units.m, 
                                                   np.sin(phi)*np.sin(theta)*I3Units.m, 
                                                   np.cos(theta)*I3Units.m)
            newPulse.time = self.flashTime
            newPulse.numberOfPhotonsNoBias = self.photonsPerPulse
            # from icecube/200704001 section 2:
            newPulse.pulseWidth = 10. * I3Units.ns
            # these two have different meanings than for flashers:
            # polar is the angle w.r.t. the candle axis,
            # azimuthal is angle along the circle (and should always be 360deg)
            newPulse.angularEmissionSigmaPolar = 180*I3Units.deg
            #newPulse.angularProfileDistributionPolar = 180*I3Units.deg
            #newPulse.angularProfileDistributionAzimuthal = 360*I3Units.deg
            newPulse.angularEmissionSigmaAzimuthal = 360.*I3Units.deg

            # insert a single pulse
            outputSeries.append(newPulse)

        frame[self.flasherPulseSeriesName] = outputSeries

        self.PushFrame(frame)
