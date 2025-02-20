'''The purpose of this library is to illustrate a way of integrating a changing detector geometry into icetray. Most functions are just placeholders without actual functionality.'''

from icecube import icetray, dataclasses, dataio
from AcousticFrame import geometry

import inspect

Acoustic = icetray.I3Frame.Stream('A')

class injectFrame(icetray.I3Module):
'''Parent module to inject new geometry, calibration and detector status frames into an existing icetray stream. To use this module, define a child that implements an inject(self) function.'''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('insertbefore',
            'Type of frame before which new frame gets injected', icetray.I3Frame.DAQ)
        self.AddParameter('delay',
            'Number of frames of type "insertbefore" before frame gets injected', 0)
        self.AddParameter('strings', 'number of strings', 1)
        self.AddParameter('modules', 'number of modules', 1)
        self.ctr = 0
    def Process(self):
        frame = self.PopFrame()
        if frame.Stop == self.GetParameter('insertbefore'):
            if self.ctr == self.GetParameter('delay'):
                self.Inject()
            self.ctr += 1
        self.PushFrame(frame)
class injectCalibrationFrame(injectFrame):
'''Module to inject new calibration frames into icetray stream.'''
    def Inject(self):
        CFrame = icetray.I3Frame(icetray.I3Frame.Calibration)
        CFrame["I3Calibration"] = dataclasses.I3Calibration()
        domcalmap   = CFrame["I3Calibration"].dom_cal
        for string in range(self.GetParameter('strings')):
            for module in range(self.GetParameter('modules')):
                omkey = icetray.OMKey(string+1, module+1)
                newdomcal = dataclasses.I3DOMCalibration()
                newdomcal.relative_dom_eff = 1.0
                domcalmap[omkey] = newdomcal
        self.PushFrame(CFrame)
class injectDetectorStatusFrame(injectFrame):
'''Module to inject new detector status frames into icetray stream.'''
    def Inject(self):
        DFrame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
        DFrame["I3DetectorStatus"] = dataclasses.I3DetectorStatus()
        domstatusmap = DFrame["I3DetectorStatus"].dom_status
        for string in range(self.GetParameter('strings')):
            for module in range(self.GetParameter('modules')):
                omkey = icetray.OMKey(string+1, module+1)
                newdomstatus = dataclasses.I3DOMStatus()
                newdomstatus.pmt_hv = 1345.*icetray.I3Units.V # arbitrary
                domstatusmap[omkey] = newdomstatus
        self.PushFrame(DFrame)
class injectAcousticFrame(injectFrame):
'''Module to inject new acoustic frames into icetray stream. These dummy acoustic frames do not contain any acoustic data, but a tilted geometry that can be copied to a geometry frame. The parameters of this module are inherited from injectFrame and createTiltGeometry.'''
    def __init__(self, context):
        injectFrame.__init__(self, context)
        source_function = geometry.createTiltedGeometry
        self.parameters = inspect.signature(source_function).parameters
        for param, sig in self.parameters.items():
            self.AddParameter(param, param, sig.default)
    def Inject(self):
        kwargs = {}
        for param, sig in self.parameters.items():
            kwargs[param] = self.GetParameter(param)
        AFrame = icetray.I3Frame(Acoustic)
        AFrame["AcousticData"] = geometry.createTiltedGeometry(**kwargs)
        self.PushFrame(AFrame)

class createGeometryFromAcoustic(icetray.I3Module):
'''Dummy module to illustrate the injection of new geometry frames after each acoustic frame.
This module just copies geometry data from the existing acoustic frames into new geometry frames.
In a more realistic scenario, this module would implement triangulation calculations based on acoustic receiver data.'''
    def Process(self):
        frame = self.PopFrame()
        self.PushFrame(frame)
        if frame.Stop == Acoustic:
            GFrame = icetray.I3Frame(icetray.I3Frame.Geometry)
            GFrame["I3Geometry"] = frame["AcousticData"]
            self.PushFrame(GFrame)

class interpolateGeometryFromAcoustic(icetray.I3Module):
'''Dummy module to illustrate geometry interpolations. Once an acoustic frame is received, all arriving frames are stored in a queue until another acoustic frame arrives. Then the queue is emptied while new geometry frames are injected. The geometry in these new geometry frames is a simple linear interpolation in x,y,x of the geometry defined in the acoustic frames.'''

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.initial_AFrame = None
        self.final_AFrame   = None
        self.frame_queue    = []
    def EmptyQueue(self):
        q_n = sum(f.Stop == icetray.I3Frame.DAQ for f in self.frame_queue)
        q_i = 0
        while self.frame_queue:
            f = self.frame_queue.pop(0)
            if f.Stop == icetray.I3Frame.DAQ:
                self.DropInterpolatedGFrame(q_i/q_n)
                q_i += 1
            self.PushFrame(f)
    def DropInterpolatedGFrame(self, scale):
        GFrame = icetray.I3Frame(icetray.I3Frame.Geometry)
        if self.initial_AFrame and self.final_AFrame:
            GFrame["I3Geometry"] = dataclasses.I3Geometry()
            geo = GFrame["I3Geometry"].omgeo
            initial_geo = self.initial_AFrame["AcousticData"].omgeo
            final_geo   = self.final_AFrame["AcousticData"].omgeo
            for key, _ in initial_geo:
                ompos = initial_geo[key].position * (1 - scale) \
                      + final_geo[key].position * scale
                newomgeo             = dataclasses.I3OMGeo()
                newomgeo.omtype      = initial_geo[key].omtype
                newomgeo.orientation = initial_geo[key].orientation
                newomgeo.position    = ompos
                geo[key]      = newomgeo
        elif self.initial_AFrame:
            GFrame["I3Geometry"] = self.initial_AFrame["AcousticData"]
        elif self.final_AFrame:
            GFrame["I3Geometry"] = self.final_AFrame["AcousticData"]
        else:
            raise('No acoustic frames available')
        self.PushFrame(GFrame)
    def Process(self):
        frame = self.PopFrame()
        if frame.Stop == Acoustic:
            self.final_AFrame = frame
            self.EmptyQueue()
            self.initial_AFrame = self.final_AFrame
            self.final_AFrame   = None
        if frame.Stop != icetray.I3Frame.Geometry:
            self.frame_queue.append(frame)
    def Finish(self):
        self.EmptyQueue()

class CreateSeparateGCD(icetray.I3Module):
'''To use this module, define your own module that inherits from this module and defines it's own InsertModule function, where the module to be wrapped is inserted into the wrapper tray.'''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        # empty dummy frames
        self.G = icetray.I3Frame(icetray.I3Frame.Geometry)
        self.C = icetray.I3Frame(icetray.I3Frame.Calibration)
        self.D = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
    def InsertModule(self, wrapper_tray, gcdfile):
        pass
    def Process(self):
        frame = self.PopFrame()
        if frame.Stop == icetray.I3Frame.Geometry:
            self.G = frame
            self.PushFrame(frame)
        elif frame.Stop == icetray.I3Frame.Calibration:
            self.C = frame
            self.PushFrame(frame)
        elif frame.Stop == icetray.I3Frame.DetectorStatus:
            self.D = frame
            self.PushFrame(frame)
        elif frame.Stop == Acoustic:
            self.PushFrame(frame)
        else:
            gcdfile = 'temporary_gcd.i3'
            outfile = dataio.I3File(gcdfile, 'w')
            outfile.push(self.G)
            outfile.push(self.C)
            outfile.push(self.D)
            outfile.close()

            class injector(icetray.I3Module):
                def __init__(self, context):
                    icetray.I3Module.__init__(self, context)
                def Process(self):
                    self.PushFrame(frame)
            outframes = []
            class extractor(icetray.I3Module):
                def __init__(self, context):
                    icetray.I3Module.__init__(self, context)
                def Process(self):
                    outframes.append(self.PopFrame())
            wrapper_tray = icetray.I3Tray()
            wrapper_tray.AddModule(injector)
            self.InsertModule(wrapper_tray, gcdfile)
            wrapper_tray.AddModule(extractor)
            wrapper_tray.Execute(1)
            wrapper_tray.Finish()
            for o in outframes:
                self.PushFrame(o)

