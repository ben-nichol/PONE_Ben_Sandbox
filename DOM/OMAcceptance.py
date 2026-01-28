import numpy as np

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey

from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.POMModel import POM



class OMAcceptance(icetray.I3ConditionalModule):
    '''
    Icetray module that separates I3Photons from CLSim into
    their respective PMTs and applies a POM acceptance cut.
    Writes the surviving photons back into the frame as I3RecoPulses
    '''

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('input_map',
                          'Name of the I3Photons from clsim',
                          'I3Photons')
        self.AddParameter('output_map',
                          'Name of output pulse series',
                          'Accepted_PulseMap')
        self.AddParameter('random_service',
                          'I3RandomService')
        self.AddParameter('drop_strings',
                          'List of string indices to ignore',
                          [])
        self.AddParameter('drop_empty',
                          'Bool to determine if empty frames should be removed from the i3 file',
                          False)


    def Configure(self):
        self.input_map      = self.GetParameter('input_map')
        self.output_map     = self.GetParameter('output_map')
        self.random_service = self.GetParameter('random_service')
        self.drop_strings   = self.GetParameter('drop_strings')
        self.drop_empty     = self.GetParameter('drop_empty')

        # load the module acceptance
        self.module = POM(random_service = self.random_service)


    def split_pmts(self, photon_map, drop_strings=[]):
        '''
        Split photon hits into PMTs on the OMs
        '''
        output_MCPE_map = simclasses.I3MCPESeriesMap() # Changed to I3MCPESeriesMap
        
        # make new map with individual PMTs
        for omkey in photon_map.keys():
            # ignore this omkey if itsplit_pmts is supposed to be dropped
            if omkey.string in drop_strings:
                continue

            new_omkey = NoPMTKey(omkey)

            # collect the photons that hit a particular OM
            #Correct as it is. new_omkey has an index for specific PMTs, something that isn't
            #assigned until the line below. omkey just cares about specific OMs
            photon_list           = np.array(photon_map[omkey])
            hit_photons, hit_pmts = self.module.apply_acceptance_cut(photon_list)

            # add these hits to the new photon map and add
            # a reco pulse to the output pulse map
            for photon, pmt in zip(hit_photons, hit_pmts):
                pmtkey = AddPMTKey(new_omkey, pmt)
                if pmtkey not in output_MCPE_map.keys():
                    output_MCPE_map[pmtkey] = simclasses.I3MCPESeriesMap()

<<<<<<< HEAD
                mcpe        = simclasses.I3MCPE()
                mcpe.time   = photon.time # Could add time delay here
                mcpe.npe    = 1 # Number of MCPEs per photon
                output_MCPE_map[pmtkey].append(mcpe)
=======
                MCPE        = dataclasses.I3MCPE()
                MCPE.time   = photon.time # Could make this non-instentaneous
                MCPE.charge = 1.0 # 1 photon = 1 PE
                output_MCPE_map[pmtkey].append(MCPE)
>>>>>>> 5ebe6e6943ae0b84b5862c415588d05627dddb8c

        return output_MCPE_map

    def DAQ(self, frame):
        photon_map = frame[self.input_map]

        # if there are no photons andsplit_pmts no noise or
        # we don't want noise just skip
        if (len(photon_map) < 1) and self.drop_empty:
            return

        frame[self.output_map] = self.split_pmts(photon_map, self.drop_strings)

        self.PushFrame(frame)
