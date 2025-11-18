#!/usr/bin/env python3
# -----------------------------
# Module to import the output of muon generator MUPAGE (ASCII table)
# for further P-ONE analysis (i3 files)
# usage: tray.AddModule(MUPAGEImporter, 'importing_muons', InFile=args.infile)
# authors: Shreya Sharma, Swathi Karanth
# -----------------------------


from icecube import icetray, dataclasses, corsika_reader
import pandas as pd


class MUPAGEImporter(icetray.I3ConditionalModule):
    """
    import the output from MUPAGE to i3 files
    """

    def __init__(self, context):
        super(MUPAGEImporter, self).__init__(context)
        self.AddParameter('InFile', 'input file name')
        self.AddParameter('OutKey', 'key of the output MCTree', 'I3MCTree')
        self.AddParameter('StoreHeader', 'stores the event header', True)

    def Configure(self):
        self.infile = self.GetParameter('InFile')
        self.outkey = self.GetParameter('OutKey')
        self.storeheader = self.GetParameter('StoreHeader')
        self.data = self.GetData(self.infile)
        self.event_nos = self.data['i_ev'].unique()
        self.event_idx = 0

    def GetData(self, filename):
        df = pd.read_csv(filename, sep=r'\s+', header=None, skiprows=3)
        df = df.rename({0: 'i_ev', 1: 'mult', 2: 'i_mu',
                        3: 'x', 4: 'y', 5: 'z',
                        6: 'dir_x', 7: 'dir_y', 8: 'dir_z',
                        9: 'energy', 10: 'time', 11: 'type'
                        }, axis=1)

        # workaround for m=1 direction bug in MUPAGE.
        # Remove once bug is solved
        df.loc[df['mult'] == 1, 'dir_z'] *= -1

        return df

    def DAQ(self, frame):
        mctree = dataclasses.I3MCTree()

        current_event = self.event_nos[self.event_idx]
        # event here refers to a bundle
        # muons inside the bundle are added as daughteres
        current_muons = self.data[self.data['i_ev']
                                  == current_event].reset_index(drop=True)

        primary = dataclasses.I3Particle()
        # populating first muon info as primary
        primary.pos.x = current_muons['x'][0]
        primary.pos.y = current_muons['y'][0]
        primary.pos.z = current_muons['z'][0]

        primary.dir = dataclasses.I3Direction(current_muons['dir_x'][0],
                                              current_muons['dir_y'][0],
                                              current_muons['dir_z'][0])

        primary.time = current_muons['time'][0]

        mctree.add_primary(primary)

        for i in range(0, len(current_muons)):
            daughter = dataclasses.I3Particle()
            daughter.pos = dataclasses.I3Position(current_muons['x'][i],
                                                  current_muons['y'][i],
                                                  current_muons['z'][i])

            daughter.dir = dataclasses.I3Direction(current_muons['dir_x'][i],
                                                   current_muons['dir_y'][i],
                                                   current_muons['dir_z'][i])

            daughter.energy = current_muons['energy'][i]
            daughter.time = current_muons['time'][i]
            # for corsika to pdg convertor use int of corsika id
            daughter.pdg_encoding = corsika_reader.CorsikaToPDG(
                int(current_muons['type'][i]))
            daughter.location_type = daughter.InIce

            mctree.append_child(primary.id, daughter)

        if self.storeheader:
            frame['I3EventHeader'] = dataclasses.I3EventHeader()
            frame['I3EventHeader'].event_id = int(current_event)

        frame[self.outkey] = mctree
        self.PushFrame(frame)
        if current_event == self.event_nos[-1]:
            self.RequestSuspension()
        self.event_idx += 1

    def Finish(self):
        print('processed events', self.event_idx)
