import numpy as np
from heapq import merge

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey

from Utilities.DOMUtility import NoPMTKey
from Utilities.PMTAcceptance import POM
from NoiseGenerators.NoiseUtility import get_mcpe_map



class DOMSimulation(icetray.I3ConditionalModule):
    '''
    Simple implementation of the POM response.

    This version uses a slightly more detailed POM
    acceptance model initially used for K40 studies
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('input_map',
                          'Name of the I3Photons from clsim',
                          'I3Photons')
        self.AddParameter('output_map',
                          'Name of the output I3RecoPulseSeriesMap',
                          'I3RecoPulseSeriesMap')
        self.AddParameter('outputmap_mcpe',
                          'Name of the output I3MCPESeriesMap',
                          'I3MCPESeriesMap')
        self.AddParameter('pmt_tts',
                          'Transit time spread of PMT',
                          3.0 * I3Units.ns)
        self.AddParameter('pmt_ts',
                          'Transit time of PMT',
                          25.0 * I3Units.ns)
        self.AddParameter('charge_sigma',
                          'Sigma of charge distribution',
                          0.3)
        self.AddParameter('charge_mean',
                          'Mean of Charge distribution',
                          1.0)
        self.AddParameter('afterpulse_prob',
                          'Total AP probability',
                          0.06)
        self.AddParameter('afterpulse_meantime_1',
                          'Mean of early time AP distribution',
                          2000.0 * I3Units.ns)
        self.AddParameter('afterpulse_timesigma_1',
                          'Sigma of early time AP distribution',
                          1000.0 * I3Units.ns)
        self.AddParameter('afterpulse_meantime_2',
                          'Mean of late time AP distribution',
                          8000.0 * I3Units.ns)
        self.AddParameter('afterpulse_timesigma_2',
                          'Sigma of late time AP distribution',
                          2000.0 * I3Units.ns)
        self.AddParameter('afterpulse_componet_ratio',
                          'Fraction of AP in early time component',
                          0.3)
        self.AddParameter('late_pulse_prob',
                          'Probability for a late pulse',
                          0.01)
        self.AddParameter('min_time_sep',
                          ' Minimum time for a separated pulse.',
                          3.0)
        self.AddParameter('pe_threshold',
                          ' Pulse charge threshold',
                          0.25)
        self.AddParameter('pe_saturation',
                          'Saturation threshold for PMT',
                          100.0)
        self.AddParameter('random_service',
                          'Random Service')
        self.AddParameter('split_doms',
                          '',
                          True)
        self.AddParameter('drop_strings',
                          '',
                          [])
        self.AddParameter('noise_pulse_series',
                          '',
                          [])
        self.AddParameter('no_pure_noise_events',
                          'Bool whether to skip events that are only noise',
                          False)
        self.AddParameter('drop_empty',
                          'Bool to determine if empty frames should be removed from the i3 file',
                          False)
        self.AddParameter('use_dark',
                          'Bool to determine if dark noise will be included',
                          False)
        self.AddParameter('use_k40',
                          'Bool to determine if  noik40se will be included',
                          False)
        self.AddParameter('dark_map',
                          'Name of the dark hit map',
                          'DarkHits')
        self.AddParameter('k40_map',
                          'Name of the k40 hit map',
                          'K40Hits')
        self.AddOutBox('OutBox')


    def Configure(self):
        self.input_map                 = self.GetParameter('input_map')
        self.output_map                = self.GetParameter('output_map')
        self.outputmap_mcpe            = self.GetParameter('outputmap_mcpe') # NOT USED!!!!!
        self.pmt_tts                   = self.GetParameter('pmt_tts')
        self.pmt_ts                    = self.GetParameter('pmt_ts')
        self.charge_sigma              = self.GetParameter('charge_sigma')
        self.charge_mean               = self.GetParameter('charge_mean')
        self.afterpulse_probability    = self.GetParameter('afterpulse_prob')
        self.afterpulse_meantime_1     = self.GetParameter('afterpulse_meantime_1')
        self.afterpulse_timesigma_1    = self.GetParameter('afterpulse_timesigma_1')
        self.afterpulse_meantime_2     = self.GetParameter('afterpulse_meantime_2')
        self.afterpulse_timesigma_2    = self.GetParameter('afterpulse_timesigma_2')
        self.afterpulse_componet_ratio = self.GetParameter('afterpulse_componet_ratio')
        self.late_pulse_probability    = self.GetParameter('late_pulse_prob')
        self.min_time_sep              = self.GetParameter('min_time_sep')
        self.pe_threshold              = self.GetParameter('pe_threshold')
        self.pe_saturation             = self.GetParameter('pe_saturation')
        self.random_service            = self.GetParameter('random_service')
        self.split_doms                = self.GetParameter('split_doms')
        self.drop_strings              = self.GetParameter('drop_strings')
        self.noise_pulse_series        = self.GetParameter('noise_pulse_series') # MAYBE CAN BE USED FOR EXTRA NOISE??? BIOLUMI (REPLACE K40 WITH THIS?)
        self.no_pure_noise_events      = self.GetParameter('no_pure_noise_events')
        self.drop_empty                = self.GetParameter('drop_empty')
        self.use_dark                  = self.GetParameter('use_dark')
        self.use_k40                   = self.GetParameter('use_k40')
        self.dark_map                  = self.GetParameter('dark_map')
        self.k40_map                   = self.GetParameter('k40_map')

        self.WroteActiveDOMsToSimFrame = False

        # load the updated module acceptance
        self.module   = POM()
        self.num_pmts = len(self.module.PMT_MATRIX)

        self.DEBUG_FRAME = 0



    def combine_ordered_lists(self, list_1, list_2):
        '''
        Combines two ALREADY ORDERED lists into a single
        ordered list preserving the total order using heapq.merge.
        Optimized version using standard library.
        
        Returns combined list and None (for backward compatibility).
        '''
        # Use heapq.merge for optimal merging of sorted iterables
        combined_list = list(merge(list_1, list_2, key=lambda x: x[0]))
        return combined_list, None


    def apply_pmt_timing_characteristics(self, mcpe_map, pmt_num, late_pulses=True, after_pulses=True):
        '''
        Applies transit time and transit time spread
        to a given mcpe map

        late_pulses and after_pulses are two booleans
        which toggle including late pulses and afterpulses
        '''
        pulse_time_list = []
        # pulseseries   = dataclasses.I3RecoPulseSeries()
        # dompulseseries = dataclasses.I3RecoPulseSeries()
        pulse_out_of_order_time_list = []

        # for every photoelectron shift the time
        # based on tts and late pulse probability
        for pe in mcpe_map:
            # time shifted by tts
            # Changed to take proper MCPE object attributes
            time = self.random_service.gaus(pe.time, self.pmt_tts)
            # check if the pulse will be a late pulse based
            # on the late pulse probability
            if late_pulses:
                if self.random_service.uniform(0.0, 1.0) < self.late_pulse_probability:
                    time += self.random_service.gaus(self.pmt_ts * 2.0, np.sqrt(2.0) * self.pmt_tts)

            if len(pulse_time_list) < 1 or time > pulse_time_list[-1][0]:
                pulse_time_list.append((time, pmt_num))
            else:
                pulse_out_of_order_time_list.append((time, pmt_num))

        # if there are pulses out of order we need to combine the
        # two ordered lists into one
        if len(pulse_out_of_order_time_list) > 0:
            pulse_out_of_order_time_list.sort(key=lambda x: x[0])

            # Use optimized heapq.merge (no need to extract times separately)
            pulse_time_list, _ = self.combine_ordered_lists(pulse_time_list,
                                                            pulse_out_of_order_time_list)

        # now add afterpulses
        if after_pulses:
            # Fixed: Initialize list BEFORE loop to collect all afterpulses
            pulse_out_of_order_time_list = []
            for pe in pulse_time_list:
                if self.random_service.uniform(0.0, 1.0) < self.afterpulse_probability:
                    if self.random_service.uniform(0.0, 1.0) < self.afterpulse_componet_ratio:
                        time = pe[0] + self.random_service.gaus(self.afterpulse_meantime_1, self.afterpulse_timesigma_1)
                    else:
                        time = pe[0] + self.random_service.gaus(self.afterpulse_meantime_2, self.afterpulse_timesigma_2)
                    pulse_out_of_order_time_list.append((time, pe[1]))

            # if there are pulses or afteruplses out of order we need
            # to combine the two ordered lists into one
            if len(pulse_out_of_order_time_list) > 0:
                pulse_out_of_order_time_list.sort(key=lambda x: x[0])

                # Use optimized heapq.merge (no need to extract times separately)
                pulse_time_list, _ = self.combine_ordered_lists(pulse_time_list,
                                                                pulse_out_of_order_time_list)

        return pulse_time_list


    def make_reco_pulse(self, pulse_time_list, pulse_charge_list, omkey, output_pulse_map, om_pulse_map=None):
        '''
        Populates the output_pulse_map with an I3RecoPulseSeries
        based on the input pulse times and charges.
        
        pulse_time_list: list of (time, pmt_num) tuples
        pulse_charge_list: list of charge values
        '''
        min_gap   = 4.0
        min_index = -1

        if len(pulse_time_list) > 100:
            leading   = 0
            following = 1
            while following < len(pulse_time_list):
                if (
                    pulse_time_list[following][0] - pulse_time_list[leading][0]
                ) < 3.0 and pulse_charge_list[leading] * pulse_charge_list[following] > 0.0:
                    pulse_charge_list[leading] += pulse_charge_list[following]
                    pulse_charge_list[following] = 0.0
                elif pulse_charge_list[following] > 0.0:
                    leading = following
                following += 1
        else:
            # needs to be better
            for i in range(1, len(pulse_time_list)):
                if (
                    pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                ) < min_gap and pulse_charge_list[i] * pulse_charge_list[i - 1] > 0.0:
                    min_gap = pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                    min_index = i
            # If less than limit, combine pulses
            while min_gap <= self.min_time_sep:
                if pulse_charge_list[min_index] > pulse_charge_list[min_index - 1]:
                    pulse_charge_list[min_index] += pulse_charge_list[min_index - 1]
                    pulse_charge_list[min_index - 1] = 0.0
                else:
                    pulse_charge_list[min_index - 1] += pulse_charge_list[min_index]
                    pulse_charge_list[min_index] = 0.0
                min_gap = self.min_time_sep + 1.0
                min_index = -1
                # reestablish new min gap
                for i in range(1, len(pulse_time_list)):
                    if (
                        pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                    ) < min_gap and pulse_charge_list[i] * pulse_charge_list[i - 1] > 0.0:
                        min_gap = pulse_time_list[i][0] - pulse_time_list[i - 1][0]
                        min_index = i

        # Optimized: Use set for O(1) lookups instead of list
        pmt_in_set = set()
        for i, (time, pmt_num) in enumerate(pulse_time_list):
            if pulse_charge_list[-1 - i] >= self.pe_threshold:
                pmt_in_set.add(pmt_num)
        pmt_in_list = sorted(pmt_in_set)

        # Simplified: Use omkey directly (no AddPMTKey needed)
        output_pulse_map[omkey] = dataclasses.I3RecoPulseSeries()

        for i in range(len(pulse_time_list)):
            # remove pulses with too low charge.
            if pulse_charge_list[-1 - i] < self.pe_threshold:
                continue
            rpulse = dataclasses.I3RecoPulse()
            rpulse.time = pulse_time_list[i][0]
            # saturate pulses with too much charge.
            if pulse_charge_list[-1 - i] > self.pe_saturation:
                rpulse.charge = self.pe_saturation
                rpulse.charge += (pulse_charge_list[-1 - i] - self.pe_saturation) * (
                    self.pe_saturation / pulse_charge_list[-1 - i]
                )
            else:
                rpulse.charge = pulse_charge_list[-1 - i]
            # The pulse width is being set to the PMT number. Don't know why that is.
            # This will need to be changed in the future to be something more physical.
            rpulse.width = pulse_time_list[i][1]
            if not (om_pulse_map is None):
                if omkey not in om_pulse_map:
                    om_pulse_map[omkey] = dataclasses.I3RecoPulseSeries()
                om_pulse_map[omkey].append(rpulse)

            output_pulse_map[omkey].append(rpulse)


    def apply_pmt_response(self, mcpe_tuple_map):
        '''
        Apply the response of the PMT, including combining pulses
        that are too close together.
        
        Now accepts mcpe_tuple_map where values are lists of (time, pmt_num) tuples.
        '''
        mcpe_tuple_map = self.apply_dead_time(mcpe_tuple_map)

        output_pulse_map = dataclasses.I3RecoPulseSeriesMap()
        om_pulse_map     = dataclasses.I3RecoPulseSeriesMap()

        # Removed .keys() - unnecessary overhead
        for omkey in mcpe_tuple_map:
            # Pre-allocate list for better performance
            mcpe_count = len(mcpe_tuple_map[omkey])
            pulse_charge_list = [self.random_service.gaus(self.charge_mean, self.charge_sigma)
                                for _ in range(mcpe_count)]
            
            self.make_reco_pulse(mcpe_tuple_map[omkey], pulse_charge_list, omkey, output_pulse_map, om_pulse_map)
        
        return output_pulse_map, om_pulse_map


    def apply_dead_time(self, mcpe_tuple_map):
        '''
        Removes successive hits on the same PMT based on the PMT dead time.
        
        Now accepts mcpe_tuple_map where values are lists of (time, pmt_num) tuples.
        Returns the same format.
        '''
        dead_time_ns          = 10.
        dead_removed_mcpe_map = {}

        for omkey in mcpe_tuple_map:
            dead_removed_mcpe_map[omkey] = []
            last_hit_times = np.ones(self.num_pmts) * -9999.

            for time_pmt_tuple in mcpe_tuple_map[omkey]:
                # Handle tuple format: (time, pmt_num)
                time = time_pmt_tuple[0]
                pmt = time_pmt_tuple[1]
                if time - last_hit_times[pmt-1] > dead_time_ns:
                    dead_removed_mcpe_map[omkey].append(time_pmt_tuple)
                    last_hit_times[pmt-1] = time
        
        return dead_removed_mcpe_map


    def collect_unique_items(self, item_lists):
        '''
        Combines all items the given lists
        and returns a list of unique items
        '''
        # Optimized: Use set for O(1) lookups instead of O(n) list lookups
        unique_set = set()
        for item_list in item_lists:
            unique_set.update(item_list)
        return list(unique_set)


    def merge_mcpe_maps(self, mcpe_map_1, mcpe_map_2):
        '''
        Combines two mcpe maps into one keeping the mcpes time ordered.
        
        Now works with tuple format where values are lists of (time, pmt_num) tuples.
        '''
        merged_map = {}

        if len(mcpe_map_1) < 1:
            return mcpe_map_2
        if len(mcpe_map_2) < 1:
            return mcpe_map_1

        # find the unique omkeys
        all_omkeys = self.collect_unique_items([list(mcpe_map_1.keys()), list(mcpe_map_2.keys())])

        for omkey in all_omkeys:
            # if the omkey is in one of the maps no need to go through the whole merge
            if omkey not in mcpe_map_1:
                merged_map[omkey] = mcpe_map_2[omkey]
                continue
            if omkey not in mcpe_map_2:
                merged_map[omkey] = mcpe_map_1[omkey]
                continue
            
            # Use optimized heapq.merge (tuples already ordered by first element)
            combined_mcpes, _ = self.combine_ordered_lists(mcpe_map_1[omkey],
                                                           mcpe_map_2[omkey])
            merged_map[omkey] = combined_mcpes
        
        return merged_map


    def Geometry(self, frame):
        # filter out the strings we want to use based
        # on dropstrings
        if len(self.drop_strings) > 0:
            self.domsUsed = dataclasses.I3OMGeoMap()
            for omkey in frame['I3Geometry'].omgeo.keys():
                if omkey.string in self.drop_strings:
                    continue
                self.domsUsed[omkey] = frame['I3Geometry'].omgeo[omkey]
        else:
            self.domsUsed = frame['I3Geometry'].omgeo


        domkeylist   = []
        self.nstring = 0
        self.nom     = 0
        self.npmt    = 0

        for omkey in self.domsUsed.keys():
            self.nstring = max(self.nstring, omkey.string)
            self.nom     = max(self.nom, omkey.om)
            self.npmt    = max(self.npmt, omkey.pmt)
            domkeylist.append(NoPMTKey(omkey))

        self.domkeys  = set(domkeylist)
        self.domkeys  = sorted(self.domkeys, key=lambda x: (x.string, x.om, x.pmt))
        self.nstring += 1
        self.nom     += 1
        self.npmt    += 1

        self.PushFrame(frame)


    def Simulation(self, frame):
        if len(self.drop_strings) > 0:
            frame['SimulatedDOMs'] = self.domsUsed

        self.WroteActiveDOMsToSimFrame = True

        self.PushFrame(frame)


    def DAQ(self, frame):
        # print(f'FRAME # {self.DEBUG_FRAME}')
        # self.DEBUG_FRAME += 1

        if not self.WroteActiveDOMsToSimFrame:
            simframe = icetray.I3Frame('S')
            self.Simulation(simframe)
        
        simulation_mcpe_map  = frame[self.input_map]

        length_noise_mcpes = 0
        if self.use_dark:
            dark_mcpe_map        = frame[self.dark_map]
            length_noise_mcpes += len(dark_mcpe_map)
        if self.use_k40:
            k40_mcpe_map         = frame[self.k40_map]
            length_noise_mcpes += len(k40_mcpe_map)
        
        # if there are no mcpes or noise
        # just drop the frame
        if self.drop_empty:
            if (len(simulation_mcpe_map) < 1) and self.no_pure_noise_events:
                return            
            if len(simulation_mcpe_map) < 1 and length_noise_mcpes < 1:
                return

        noise_mcpe_maps = []

        # OPTIMIZED: Keep data as tuples throughout - no I3MCPE object recreation
        # Convert I3MCPESeries to tuple format once, process as tuples
        simulation_tuple_map = {}
        for omkey in simulation_mcpe_map:
            simulation_tuple_map[omkey] = self.apply_pmt_timing_characteristics(
                simulation_mcpe_map[omkey], omkey.pmt)
        
        if self.use_dark:
            dark_tuple_map = {}
            for omkey in dark_mcpe_map:
                dark_tuple_map[omkey] = self.apply_pmt_timing_characteristics(
                    dark_mcpe_map[omkey], omkey.pmt, late_pulses=False)
            noise_mcpe_maps.append(dark_tuple_map)
        
        if self.use_k40:
            k40_tuple_map = {}
            for omkey in k40_mcpe_map:
                k40_tuple_map[omkey] = self.apply_pmt_timing_characteristics(
                    k40_mcpe_map[omkey], omkey.pmt, late_pulses=False)
            noise_mcpe_maps.append(k40_tuple_map)

        noise_tuple_map = {}
        # merge all the noise
        for i in range(len(noise_mcpe_maps)):
            noise_tuple_map = self.merge_mcpe_maps(noise_tuple_map, noise_mcpe_maps[i])

        # merge the signal and noise
        total_tuple_map = self.merge_mcpe_maps(simulation_tuple_map, noise_tuple_map)

        # apply the PMT pulse response without any noise
        (no_noise_output_pulses, _) = self.apply_pmt_response(simulation_tuple_map)

        # apply the PMT pulse response to the combined signal and noise
        (output_pulses, om_pulses) = self.apply_pmt_response(total_tuple_map)

        frame[self.output_map]              = output_pulses
        frame[self.output_map + '_nonoise'] = no_noise_output_pulses
        frame['triggerpulsemap']           = om_pulses

        self.PushFrame(frame)