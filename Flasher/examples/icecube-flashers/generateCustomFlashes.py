from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

from I3Tray import I3Units

class CustomFlashes(icetray.I3ConditionalModule):
    """
    Generates I3CLSimFlasherPulse objects for 370nm and 405nm LEDs (I&II)

    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3CLSimFlasherPulseSeries to write",
                          "I3CLSimFlasherPulseSeries")
        self.AddParameter("PhotonsPerPulse",
                          "Photons to emit per pulse (see http://wiki.icecube.wisc.edu/index.php/Standard_Candle\n" +
                          "for some more information on what number to use)",
                          5)
        self.AddParameter("FlashTime",
                          "Time (within each event) at which to flash",
                          0.*I3Units.ns)
        self.AddParameter("CandleNumber",
                          "Choose which source to simulate: 1 for LED370nm or 2 for LED405nm",
                          1)
        #self.AddParameter("NumPulses",
        #"Select the number of pulses to be generated",
        #                  1)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.flasherPulseSeriesName = self.GetParameter("FlasherPulseSeriesName")
        self.photonsPerPulse = self.GetParameter("PhotonsPerPulse")
        self.flashTime = self.GetParameter("FlashTime")
        self.candleNumber = self.GetParameter("CandleNumber")
        #self.numPulses = self.GetParameter("NumPulses")

        if self.candleNumber not in [1,2]:
            raise RuntimeError("You have to select either LED 1 or 2. You chose %u." % self.candleNumber)

    def DAQ(self, frame):
        outputSeries = I3CLSimFlasherPulseSeries()

        numPhotons = self.photonsPerPulse

        newPulse = I3CLSimFlasherPulse()
        if self.candleNumber==1:
            newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED370nm
            # from http://wiki.icecube.wisc.edu/index.php/File:Sc_geometry.jpg
            #newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 107.66*I3Units.m)
            newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 100*I3Units.m)
            newPulse.dir = dataclasses.I3Direction(0.,0., -1.) # facing down
        elif self.candleNumber==2:
            newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
            #newPulse.wavelength = 405 * I3Units.nanometer
            # from http://wiki.icecube.wisc.edu/index.php/Standard_Candle#Standard_Candle_Mark_II
            #newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 107.66*I3Units.m)
            newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 100*I3Units.m)
            newPulse.dir = dataclasses.I3Direction(0., 0., -1) # facing down
        else:
            raise RuntimeError("invalid candle number. logic error.")

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
