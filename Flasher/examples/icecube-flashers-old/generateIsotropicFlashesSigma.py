from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

from I3Tray import I3Units

import numpy as np

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
        self.AddParameter("CandlePosition",
                           "Choose wher the source needs to be placed in Z",
                           107.66)
        #self.AddParameter("NumPulses",
        #"Select the number of pulses to be generated",
        #                  1)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.flasherPulseSeriesName = self.GetParameter("FlasherPulseSeriesName")
        self.photonsPerPulse = self.GetParameter("PhotonsPerPulse")
        self.flashTime = self.GetParameter("FlashTime")
        self.candleNumber = self.GetParameter("CandleNumber")
        self.candlePosition = self.GetParameter("CandlePosition")
        #self.numPulses = self.GetParameter("NumPulses")

        if self.candleNumber not in [1,2,3]:
            raise RuntimeError("You have to select either LED 1 or 2. You chose %u." % self.candleNumber)

    def DAQ(self, frame):
        outputSeries = I3CLSimFlasherPulseSeries()
        numPhotons = self.photonsPerPulse
        #print outputSeries
        #azi = np.arange(0, 359, 1)
        azi = np.array([0])
        czen = np.linspace(1, -1, 500000)
        print np.degrees(np.arccos(czen))
        print np.degrees(np.arccos(czen[1:]) - np.arccos(czen[:-1]))
        newPulse = ([])
        index = 0

        for z in range(0, len(czen)):
            #print(z)
            for a in range(0, len(azi)):
                #print(a)
                newPulse = np.append(newPulse, I3CLSimFlasherPulse())
                #print index
                newPulse[index].pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, self.candlePosition*I3Units.m)
                #newPulse[index].pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 107.66*I3Units.m)
                #newPulse[index].pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 100.0*I3Units.m)

                cosDiff = abs(czen[1] - czen[0])/2 - 1e-8
                if z == 0:
                    #continue
                    thetaDiff = np.arccos(czen[z]-cosDiff)
                    #thetaDiff2 = np.arccos(-1) - thetaDiff
                    newPulse[index].angularEmissionSigmaPolar = thetaDiff*I3Units.rad
                    #print thetaDiff*57
                    #print czen[z]-cosDiff

                if z == len(czen) - 1:
                    #continue
                    thetaDiff = np.arccos(czen[z]+cosDiff)
                    thetaDiff2 = np.arccos(-1) - thetaDiff
                    newPulse[index].angularEmissionSigmaPolar = thetaDiff2*I3Units.rad
                    #print thetaDiff2*57
                    #print czen[z]+cosDiff

                if z != 0 and z != len(czen) - 1:
                    #else:
                    thetaDiff = abs(np.arccos(czen[z]+cosDiff) - np.arccos(czen[z]-cosDiff))
                    newPulse[index].angularEmissionSigmaPolar = thetaDiff*I3Units.rad

                newPulse[index].dir = dataclasses.I3Direction(np.arccos(czen[z])*I3Units.rad, azi[a]*I3Units.deg) # facing down
                if self.candleNumber==1:
                    newPulse[index].type = I3CLSimFlasherPulse.FlasherPulseType.LED370nm
                    # from http://wiki.icecube.wisc.edu/index.php/File:Sc_geometry.jpg
                    #newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 107.66*I3Units.m)
                elif self.candleNumber==2:
                    newPulse[index].type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
                    #newPulse.wavelength = 405 * I3Units.nanometer
                    # from http://wiki.icecube.wisc.edu/index.php/Standard_Candle#Standard_Candle_Mark_II
                    #newPulse.pos = dataclasses.I3Position(0*I3Units.m, 0*I3Units.m, 107.66*I3Units.m)
                elif self.candleNumber==3:
                    newPulse[index].type = I3CLSimFlasherPulse.FlasherPulseType.LED450nm
                else:
                    raise RuntimeError("invalid candle number. logic error.")

                newPulse[index].time = self.flashTime
                newPulse[index].numberOfPhotonsNoBias = self.photonsPerPulse

                # from icecube/200704001 section 2:
                newPulse[index].pulseWidth = 10. * I3Units.ns

                # these two have different meanings than for flashers:
                # polar is the angle w.r.t. the candle axis,
                # azimuthal is angle along the circle (and should always be 360deg)
                #newPulse[index].angularEmissionSigmaPolar = 0.4*I3Units.deg
                #newPulse[index].angularEmissionSigmaAzimuthal = 0*I3Units.deg
                newPulse[index].angularEmissionSigmaAzimuthal = 360.*I3Units.deg

                # insert multiple pulses
                outputSeries.append(newPulse[index])
                index = index + 1

        #print outputSeries
        frame[self.flasherPulseSeriesName] = outputSeries

        self.PushFrame(frame)
