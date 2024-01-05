from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
from scipy import stats
from iminuit import minimize
from scipy.stats.distributions import chi2
from Utility.NuTauPeakShapes import log_likelihood_biGauss, log_likelihood_doublePeak
from Utility.NuTauPeakShapes import likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
from Utility.NuTauPeakShapes import log_likelihood_expGauss, log_likelihood_expDoublePeak, expGauss, expDoublePeak
import scipy, csv
#from tabulate import tabulate

class curveFit(icetray.I3ConditionalModule):
    """
    Fitting single and double bifurcated gaussian
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("omgeo",
                        "geometry map given for the analysis",
                        "omgeo")

        self.AddParameter("InputMCPETree",
                         "Input MCPETree name for analysis",
                         "MCPESeriesMap")

        self.AddParameter("OutputMCPETree",
                         "Output MCPETree name",
                         "CurvefitParameters")

        self.AddParameter("HitsInDOMsCut",
                          "Cut in the num fits in DOMS",
                          200)

        self.AddParameter("FrameList",
                          "List of frame numbers to debug",
                          [])

        self.AddParameter("StringList",
                          "List of string numbers to debug",
                          [])

        self.AddParameter("DOMList",
                          "List of DOM numbers to debug",
                          [])

        self.AddParameter("waveformmap","Name of the Waveform map series","")

        self.AddOutBox("OutBox")

    def Configure(self):

        self.omgeo = self.GetParameter("omgeo")
        self.input = self.GetParameter("InputMCPETree")
        self.output = self.GetParameter("OutputMCPETree")
        self.cuts = self.GetParameter("HitsInDOMsCut")
        self.frames = self.GetParameter("FrameList")
        self.strings = self.GetParameter("StringList")
        self.doms = self.GetParameter("DOMList")
        self.waveformmapname = self.GetParameter("waveformmap")

    def Geometry(self,frame):
        self.nopmtkeys = set()
        self.domsUsed = frame['I3Geometry'].omgeo
        self.nstring = 0
        self.nom = 0
        self.npmt = 0
        for omkey in self.domsUsed.keys() :
            self.nstring = max(self.nstring,omkey.string)
            self.nom = max(self.nom,omkey.om)
            self.npmt = max(self.npmt,omkey.pmt)
            self.nopmtkeys.add(NoPMTKey(omkey))
        self.nstring += 1
        self.nom += 1
        self.npmt += 1
        self.PushFrame(frame)


    def DAQ(self, frame):

        recoPulseMap = frame[self.input]
        waveformmapseries = frame[self.waveformmapname] 

        biGauss_valuesMap = dataclasses.I3MapKeyVectorDouble()
        doublePeak_valuesMap = dataclasses.I3MapKeyVectorDouble()
        exitStatusMap = dataclasses.I3MapKeyVectorDouble()

        for omkey in self.nopmtkeys :
            waveform = dataclasses.I3Waveform()
            mintime = 10000000.
            maxtime = 0.
            binwidth = 4.0
            for ipmt in range(self.npmt):
                omkey_wpmt = AddPMT(omkey,ipmt)
                if omkey_wpmt in waveformmapseries.keys() :
                    for wave in waveformmapseries[omkey_wpmt]:
                        mintime = min(mintime,wave.time)
                        maxtime = max(maxtime,wave.time+wave.binwidth*len(wave.waveform))
                        binwidth = wave.binwidth
            wavearray = []
            for i in range((maxtime-mintime)/binwidth):
                wavearray.append(0.0)
            for ipmt in range(self.npmt):
                omkey_wpmt = AddPMT(omkey,ipmt)
                if omkey_wpmt in waveformmapseries.keys() :
                    for wave in waveformmapseries[omkey_wpmt]:
                        startbin = (wave.time-mintime)/binwidth
                        for i in range(len(wave.waveform)):
                            wavearray[startbin+i] += wave.waveform[i]

            totalcharge = sum(wavearray)
            if totalcharge > -100.:
                exit_status = np.array([0])
                exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})

            '''
            Calculating the mean and removing the tails
            '''
            
            mean = sum([wavearray[i]*(i+1)/totalcharge for i in range(len(wavearray))])-1
            mean_charge = sum([wavearray[i] for i in range(max(0,int(mean-50./binwidth)),max(len(wavearray),int(mean+50./binwidth)))])

            if len(mean_charge) > -20:
                exit_status = np.array([1])
                exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})
                continue

            #Shifting mean to zero
            max_hitTimes_mean = sum(max_hitTimes*max_charge)/sum(max_charge)
            timestamps = max_hitTimes - max_hitTimes_mean
            final_mean = sum(timestamps*max_charge)/sum(max_charge)

            '''
            Histogramming the data from simulation
            '''
            bins = np.arange(min(timestamps), max(timestamps), 3)
            num, bin_edges = np.histogram(timestamps, bins=bins, weights=max_charge)
            bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

            num_ampRatio = num/max(num)

            #removing bins which are <1/5 the max(num), removing the tails this way.
            num_select = num[num_ampRatio > 0.0]
            bin_centers_select = bin_centers[num_ampRatio > 0.0]

            '''
            Including continuity in the bins
            '''

            #considering two extra bins on both sides
            bin_center_bool = (bin_centers >= min(bin_centers_select) - 6)&(bin_centers <= max(bin_centers_select) + 6)
            entries_in_bins = num[bin_center_bool]
            bin_centers = bin_centers[bin_center_bool]

            '''
            Removing DOMs which don't have enough hits
            '''

            if len(entries_in_bins) < 9:
                exit_status = np.array([2])
                exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})
                continue


            if max(bin_centers) <= 0:
                maxBinCenter = max(bin_centers) + abs(max(bin_centers)) + 3
            else:
                maxBinCenter = max(bin_centers)

            exitStatusMap.update({omkey: dataclasses.I3VectorDouble(exit_status)})
            time_window = max(bin_centers) - min(bin_centers)
            exit_status = np.array([3])

            '''
            Fitting bifurcated Gaussian and double bifurcated gaussian to
            the mcpe hit time distributions for both tau and electron.

            initial_biGauss = np.array([final_mean, time_window/2, 1, max(entries_in_bins)])
            initial_doublePeak = np.array([peak_time_boundary-1,
                                           20,
                                           1,
                                           max(entries_in_bins),
                                           peak_time_boundary+1,
                                           20,
                                           1,
                                           max(entries_in_bins)])


            '''

            init_k = np.linspace(1, 20, 1)
            init_wid = np.linspace(0, time_window, 20)

            best_fcn_single = 1e12
            soln_biGauss = 0
            bnds_biGauss = [[min(bin_centers), maxBinCenter],
                            [1, time_window], # Let the width be negative
                            [1, 20], # Restrict k to be positive, but only up to 20
                            [0.1, 2*max(entries_in_bins)]] # Don't restrict the amplitude, it will vary greatly with K

            for iK in init_k:
                for iwid in init_wid:
                    initial_biGauss = np.array([final_mean, iwid, iK, max(entries_in_bins)])

                    #Single Peak
                    nll = lambda *args: log_likelihood_biGauss(*args)

                    soln_single = minimize(log_likelihood_expGauss, initial_biGauss,
                                            args=(entries_in_bins, bin_centers, debug_mode),
                                            bounds = bnds_biGauss)

                    if soln_single.fun < best_fcn_single:
                        best_fcn_single = soln_single.fun
                        soln_biGauss = soln_single

            #Double Peak
            peak_time_boundary = final_mean-6.
            bnds_doublePeak = [[min(bin_centers), peak_time_boundary],
                                bnds_biGauss[1],
                                bnds_biGauss[2],
                                bnds_biGauss[3],
                                [peak_time_boundary, maxBinCenter],
                                bnds_biGauss[1],
                                bnds_biGauss[2],
                                bnds_biGauss[3]]

            best_fcn_double = 1e12
            soln_doublePeak = 0

            for iK in init_k:
                for iwid in init_wid:
                    initial_doublePeak = np.array([peak_time_boundary-1,
                                                   iwid,
                                                   iK,
                                                   max(entries_in_bins),
                                                   peak_time_boundary+1,
                                                   iwid,
                                                   iK,
                                                   max(entries_in_bins)])

                    nll = lambda *args: log_likelihood_doublePeak(*args)
                    if debug_mode == True:

                    	soln_double = minimize(log_likelihood_expDoublePeak, initial_doublePeak,
                                                args=(entries_in_bins, bin_centers, debug_mode),
                                                bounds=bnds_doublePeak)

                    if soln_double.fun < best_fcn_double:
                        best_fcn_double = soln_double.fun
                        soln_doublePeak = soln_double

            '''
            Calculating the Likelihood ratio for bifurcated gaussian
            and double double bifurcated gaussian
            '''
            LR_biGauss = likelihood_ratio_biGauss(bin_centers, entries_in_bins, soln_biGauss.x[0],
                                                  soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
            LR_doublePeak = likelihood_ratio_doublePeak(bin_centers, entries_in_bins, soln_doublePeak.x[0],
                                                        soln_doublePeak.x[1],soln_doublePeak.x[2],
                                                        soln_doublePeak.x[3], soln_doublePeak.x[4],
                                                        soln_doublePeak.x[5], soln_doublePeak.x[6],
                                                        soln_doublePeak.x[7])
            '''
            Update values
            '''

            biGauss_values = np.array([soln_biGauss.fun, soln_biGauss.x[0],
                                            soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3],
                                            LR_biGauss])

            doublePeak_values = np.array([soln_doublePeak.fun, soln_doublePeak.x[0],
                                        soln_doublePeak.x[1],soln_doublePeak.x[2],
                                        soln_doublePeak.x[3], soln_doublePeak.x[4],
                                        soln_doublePeak.x[5], soln_doublePeak.x[6],
                                        soln_doublePeak.x[7], LR_doublePeak])

            # Define omkey:vector dictionary
            biGauss_valuesMap.update({omkey: dataclasses.I3VectorDouble(biGauss_values)})
            doublePeak_valuesMap.update({omkey: dataclasses.I3VectorDouble(doublePeak_values)})

        frame[self.output+'_biGauss'] = biGauss_valuesMap
        frame[self.output+ '_doublePeak'] = doublePeak_valuesMap
        frame[self.output+ '_exitStatus'] = exitStatusMap
        # Increase the frame counter
        self.frame_counter += 1

        self.PushFrame(frame)
