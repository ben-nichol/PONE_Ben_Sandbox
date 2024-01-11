#!/usr/bin/env python
# This is meant to be a slightly more robust approach to reconstruction of a muon event.
# The physics and likelihood model is heavily based off of the ICECUBE model and can be found at
# "https://publications.ub.uni-mainz.de/theses/volltexte/2014/3869/pdf/3869.pdf"
# The time residuals are computed by myself though. The techniques used are detailed in a text document I have somewhere -dg

# Import some useful ICECUBE modules
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
from icecube.dataclasses import I3Particle
from icecube.phys_services import I3Calculator as calc
import numpy as np
import time
from Utilities.PandelPDFs import cpandel as pdf  # This module is used to store the pdf
from scipy import optimize as op
import sys
import argparse
import math as m
from Utilities.RecoUtility import GetGeoTime
from Utilities.DOMUtility import NoPMTKey, AddPMTKey
from Utilities.OpticalParameters import c, n, ngroup, tau


# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(data, domsUsed, track):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    dir = track.dir
    vertex = track.pos

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(darkAmp, StartL, StopL):
        sum_nloglike = 0.0
        t0 = track.time
        for dom in pulse_series.keys():
            em_pos = calc.cherenkov_position(track, geo_doms[dom].position, ngroup, n)
            vertex_to_em = dataclasses.I3Position(
                em_pos.x - vertex.x, em_pos.y - vertex.y, em_pos.z - vertex.z
            )
            sign = (
                vertex_to_em.x * dir.x + vertex_to_em.y * dir.y + vertex_to_em.z * dir.z
            ) / np.sqrt(
                vertex_to_em.x**2.0 + vertex_to_em.y**2.0 + vertex_to_em.z**2.0
            )
            Length_Along_Track = sign * np.sqrt(
                vertex_to_em.x**2.0 + vertex_to_em.y**2.0 + vertex_to_em.z**2.0
            )
            pos = geo_doms[dom].position
            _dir = track.dir
            vert = vertex
            d, dc, t = GetGeoTime(
                [pos.x, pos.y, pos.z],
                [vert.x, vert.y, vert.z],
                [_dir.x, _dir.y, _dir.z],
            )
            p_charge = np.exp(-d / tau) / max(dc, 0.5)
            for pulse in pulse_series[dom]:
                time_r = pulse.time - t0 - t
                cpandel_out = pdf(time_r, d)
                charge = pulse.charge
                if Length_Along_Track < StartL or Length_Along_Track > StopL:
                    sum_nloglike -= charge * np.log(darkAmp)
                else:
                    sum_nloglike -= charge * np.log(
                        (1.0 - darkAmp) * cpandel_out * p_charge + darkAmp
                    )

        return sum_nloglike

    return likelihoodFunction


class StartStopFit(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter(
            "pulseseries", "Name of the Merged MCPE tree name", "MergedSeriesMap"
        )
        self.AddParameter("seedtrack", "Track to seed fit", "linefit")
        self.AddParameter("output", "Track to store fit.", "llnfit")
        self.AddOutBox("OutBox")

    def Configure(self):
        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.MaxDOMRad = 0.0
        self.MaxZ = 0.0

        self.domsUsed = {}

    def Geometry(self, frame):
        self.domsUsed = frame["I3Geometry"].omgeo

        for dom in self.domsUsed.keys():
            dompos = self.domsUsed[dom].position
            radius = np.sqrt(dompos.x**2.0 + dompos.y**2.0)
            if radius > self.MaxDOMRad:
                self.MaxDOMRad = radius
            if abs(dompos.z) > self.MaxZ:
                self.MaxZ = dompos.z

        self.PushFrame(frame)

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.

    def Physics(self, frame):
        if not frame.Has(self.pulseseries):
            self.PushFrame(frame)
            return
        data = frame[self.pulseseries]

        if not frame.Has(self.seedtrack):
            self.PushFrame(frame)
            return
        linefit = frame[self.seedtrack]

        qFunctor = LikelihoodFunctor(data, self.domsUsed, linefit)

        # Minimize using scipy
        def func(x):
            darkAmp, startL, stopL = x
            return qFunctor(darkAmp, startL, stopL)

        solution = op.minimize(
            fun=func, x0=np.array([0.5, 0.0, 3000.0]), method="Nelder-Mead"
        )

        darkAmp = solution.x[0]
        startL = max(0.0, solution.x[1])
        stopL = solution.x[2]
        stop_pos_x = linefit.pos.x + stopL * linefit.dir.x
        stop_pos_y = linefit.pos.y + stopL * linefit.dir.y
        stop_pos_z = linefit.pos.z + stopL * linefit.dir.z
        stop_cylr = np.sqrt(stop_pos_x**2.0 + stop_pos_y**2.0)

        if abs(stop_pos_z) > self.MaxZ or stop_cylr > self.MaxDOMRad:
            l = stopL
            if abs(stop_pos_z) > 0.0 and abs(linefit.dir.z) > 0.0:
                l = (
                    (stop_pos_z / abs(stop_pos_z)) * self.MaxZ - linefit.pos.z
                ) / linefit.dir.z
            lrad = stopL
            dir_rad = np.sqrt(linefit.dir.x**2.0 + linefit.dir.y**2.0)

            if dir_rad > 0.0 and stop_cylr > self.MaxDOMRad:
                lrad = -(linefit.pos.x * linefit.dir.x)
                lrad -= linefit.pos.y * linefit.dir.y
                sqrt_inside = max(
                    0.0,
                    lrad**2.0
                    - 4.0
                    * (
                        linefit.pos.x**2.0
                        + linefit.pos.x**2.0
                        - self.MaxDOMRad**2.0
                    )
                    * (dir_rad**2.0),
                )
                lrad += np.sqrt(sqrt_inside)
                lrad /= dir_rad**2.0
            stopL = min(l, lrad)

        # include both linefit and improved recos for comparison
        frame[self.output + "_darkAmp"] = dataclasses.I3Double(darkAmp)
        frame[self.output + "_startL"] = dataclasses.I3Double(startL)
        frame[self.output + "_stopL"] = dataclasses.I3Double(stopL)

        self.PushFrame(frame)
