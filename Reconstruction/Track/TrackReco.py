#!/usr/bin/env python
# This is meant to be a slightly more robust approach to reconstruction of a muon event.
# The physics and likelihood model is heavily based off of the ICECUBE model and can be found at
# "https://publications.ub.uni-mainz.de/theses/volltexte/2014/3869/pdf/3869.pdf"
# The time residuals are computed by myself though. The techniques used are detailed in a text document I have somewhere -dg

# Import some useful ICECUBE modules
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
from icecube.dataclasses import I3Particle
import numpy as np
import time
from Utilities.PandelPDFs import cpandel as pdf  # This module is used to store the pdf
from scipy import optimize as op
import sys, os
import argparse
import math as m
from Utilities.RecoUtility import GetGeoTime
from Utilities.DOMUtility import NoPMTKey, AddPMTKey, DOMProperties
from Utilities.OpticalParameters import c, n, ngroup, tau
import pickle
from scipy.interpolate import UnivariateSpline


# Functional that is fed data from InitialGuess for PMT locations and the PDF we wish to use. Uses those locations to build a Pandel Function for a given track
def LikelihoodFunctor(_sigmascale, data, domsUsed, _vertexrad, _tau):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    sigmascale = _sigmascale
    vertexRad = _vertexrad
    pdark = 1e-14
    # min time index for the first hit PMT
    sigma = 1.40291905
    lambda_s = 2.02325418e01
    rho = 1.47088227e-02
    tau = _tau
    inputdict = pickle.load(
        open(os.getenv("PONESRCDIR") + "/data/PandelNormalization.pkl", "rb")
    )
    normalize = UnivariateSpline(inputdict["dist"], inputdict["normalization"], k=2)

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(x):
        vtheta, vphi, theta, phi, t0, Energy = x
        vertex = [
            vertexRad * np.sin(vtheta) * np.cos(vphi),
            vertexRad * np.sin(vtheta) * np.sin(vphi),
            vertexRad * np.cos(vtheta),
        ]
        direction = [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
        sum_nloglike = 0.0
        for dom in pulse_series.keys():
            d, dc, t, theta, phi = GetGeoTime(
                [
                    geo_doms[dom].position.x,
                    geo_doms[dom].position.y,
                    geo_doms[dom].position.z,
                ],
                vertex,
                direction,
            )
            pcharge = np.exp(-d / tau) / min(1.0, dc)
            for pulse in pulse_series[dom]:
                charge = 1.0
                time_r = pulse.time - t0 - t
                cpandel_out = pdf(
                    time_r,
                    d,
                    np.sqrt(sigma * sigma + sigmascale * sigmascale),
                    lambda_s,
                    rho,
                ) / normalize(min(max(1.0, d), 200.0))
                if type(pulse_series) == "icecube.dataclasses.I3RecoPulseSeriesMap":
                    charge = pulse.charge
                # print("dist = "+str(d)+" time_r = "+str(time_r)+" cpandel_out = "+str(cpandel_out))
                sum_nloglike -= charge * (np.log(Energy * cpandel_out + pdark))
                # if time_r < 0 :
                #    sum_nloglike -= charge*(np.log(Energy*cpandel_out*pcharge+pdark) + time_r)
                # else :
                #    sum_nloglike -= charge*np.log(Energy*cpandel_out*pcharge+pdark)

        return sum_nloglike

    return likelihoodFunction


def GetVertexTime(vertex, direction, pulse_series, geo_doms):
    totalcharge = 0.0
    MaxChargeDOM = None
    maxCharge = 0.0
    DOMCharge = {}

    c_n = c / ngroup  # light in water
    theta_c = np.arccos(1.0 / n)
    for dom in pulse_series.keys():
        totalcharge = 0.0
        for pulse in pulse_series[dom]:
            totalcharge += pulse.charge
        if NoPMTKey(dom) in DOMCharge.keys():
            DOMCharge[NoPMTKey(dom)] += totalcharge
        else:
            DOMCharge[NoPMTKey(dom)] = totalcharge

    # time of largest pulse
    maxCharge = 0.0
    maxCharge_time = 0.0
    for dom in DOMCharge.keys():
        if DOMCharge[dom] > maxCharge:
            maxCharge = DOMCharge[dom]
            MaxChargeDOM = dom

    if type(MaxChargeDOM) != type(OMKey(0, 0, 0)):
        return 7200

    DOMPos = geo_doms[AddPMTKey(MaxChargeDOM, 1)].position

    maxCharge = 0.0
    maxCharge_time = 0.0
    for domkey in pulse_series.keys():
        if (domkey.string == MaxChargeDOM.string) and (domkey.om == MaxChargeDOM.om):
            for pulse in pulse_series[domkey]:
                if pulse.charge > maxCharge:
                    maxCharge = pulse.charge
                    maxCharge_time = pulse.time

    x = DOMPos.x - vertex.x
    y = DOMPos.y - vertex.y
    z = DOMPos.z - vertex.z
    # Compute (\vec{r} - vec{x}) dot \vec{v}
    dotprod = x * direction.x + y * direction.y + z * direction.z
    # Compute the final vector components
    # Compute t_i,c and d_i,c
    sqrt_inside = max(0.0, x * x + y * y + z * z - dotprod * dotprod)
    dc = np.sqrt(sqrt_inside)
    # time to travel to closest approach
    tc = dotprod / c

    # Now we find the time of the photon emission
    _tc = tc - dc / (np.tan(theta_c) * c)
    # The first component of the geometric time
    d = dc / np.sin(theta_c)
    t_geo = d / c_n
    # The total geometric time
    t_geo = t_geo + _tc
    # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
    return maxCharge_time - t_geo


class TrackReco(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter(
            "pulseseries", "Name of the Merged MCPE tree name", "MergedSeriesMap"
        )
        self.AddParameter("seedtrack", "Track to seed fit", ["linefit"])
        self.AddParameter("output", "Track to store fit.", "llnfit")
        self.AddParameter("vertexRad", "Radius to put vertex at", 200.0)
        self.AddParameter("UseMC", "Use MC Truth Track to seed", False)
        self.AddParameter("tau", "Optical attenuation length", tau)
        self.AddParameter("minr", "Minimum radius for attenuation", 0.25)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.pulseseries = self.GetParameter("pulseseries")
        self.seedtrack = self.GetParameter("seedtrack")
        self.output = self.GetParameter("output")
        self.vertexRad = self.GetParameter("vertexRad")
        self.useMC = self.GetParameter("UseMC")
        self.domsUsed = {}
        # Some quantities that are environment dependent
        self.theta_c = np.arccos(1.0 / n)  # Cherenkov angle in water in radians
        self.tau = self.GetParameter("tau")
        self.minr = self.GetParameter("minr")
        self.domprop = DOMProperties()
        self.data = None

    # Main function of this file. Structured this way so that it can be easily imported aswell in any other implementation.

    def Geometry(self, frame):
        self.domsUsed = frame["I3Geometry"].omgeo

        maxradius = 0.0
        for dom in self.domsUsed.keys():
            pos = self.domsUsed[dom].position
            radius = np.sqrt(pos.x**2.0 + pos.y**2.0 + pos.z**2.0)
            maxradius = max(maxradius, radius)

        self.vertexRad += maxradius

        self.PushFrame(frame)

    def Physics(self, frame):
        if not frame.Has(self.pulseseries):
            self.PushFrame(frame)
            return

        data = frame[self.pulseseries]

        results = []
        result_nlogl = []

        qFunctor_15 = LikelihoodFunctor(15.0, data, self.domsUsed, self.vertexRad, 50.0)
        qFunctor_10 = LikelihoodFunctor(10.0, data, self.domsUsed, self.vertexRad, 40.0)
        qFunctor_5 = LikelihoodFunctor(5.0, data, self.domsUsed, self.vertexRad, 30.0)
        qFunctor_0 = LikelihoodFunctor(0.0, data, self.domsUsed, self.vertexRad, 28.0)

        functionlist = [qFunctor_15, qFunctor_10, qFunctor_5, qFunctor_0]

        for seedtrack in self.seedtrack:
            if not frame.Has(seedtrack):
                continue

            linefit = frame[seedtrack]

            direction = dataclasses.I3Direction(
                linefit.dir.x, linefit.dir.y, linefit.dir.z
            )
            p_2 = linefit.pos.x**2.0 + linefit.pos.y**2.0 + linefit.pos.z**2.0
            pd = (
                linefit.pos.x * direction.x
                + linefit.pos.y * direction.y
                + linefit.pos.z * direction.z
            )
            r_2 = self.vertexRad**2.0

            if p_2 > self.vertexRad:
                ratio = self.vertexRad / p_2
                vertex = dataclasses.I3Position(
                    linefit.pos.x * ratio, linefit.pos.y * ratio, linefit.pos.z * ratio
                )
            elif pd**2.0 - p_2 + r_2 > 0.0:
                L = -pd - np.sqrt(pd**2.0 - p_2 + r_2)
                vertex = dataclasses.I3Position(
                    linefit.pos.x + L * direction.x,
                    linefit.pos.y + L * direction.y,
                    linefit.pos.z + L * direction.z,
                )

            VTheta = vertex.theta
            VPhi = vertex.phi

            T0 = linefit.time
            if self.seedtrack == "linefit":
                T0 = GetVertexTime(vertex, direction, data, self.domsUsed)

            solution = None
            seedvalues = [VTheta, VPhi, direction.theta, direction.phi, T0, 1.0]
            for fitfunc in functionlist:
                solution = op.minimize(
                    fun=fitfunc, x0=np.array([seedvalues]), method="Nelder-Mead"
                )
                seedvalues = solution.x

            results.append(solution)
            result_nlogl.append(solution.fun)

        min_index = result_nlogl.index(min(result_nlogl))
        vx = (
            self.vertexRad
            * np.sin(results[min_index].x[0])
            * np.cos(results[min_index].x[1])
        )
        vy = (
            self.vertexRad
            * np.sin(results[min_index].x[0])
            * np.sin(results[min_index].x[1])
        )
        vz = self.vertexRad * np.cos(results[min_index].x[0])
        q = dataclasses.I3Position(vx, vy, vz)
        phi = results[min_index].x[3]
        theta = results[min_index].x[2]
        u = dataclasses.I3Direction(
            np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
        )

        # Record the final result
        recoParticle = dataclasses.I3Particle()
        recoParticle.shape = dataclasses.I3Particle.InfiniteTrack

        # record on particle whether reconstruction was successful
        if results[min_index].success == True:
            recoParticle.fit_status = dataclasses.I3Particle.OK
        else:
            recoParticle.fit_status = dataclasses.I3Particle.InsufficientQuality

        recoParticle.dir = u
        recoParticle.speed = c
        recoParticle.pos = q
        recoParticle.time = results[min_index].x[4]
        # print("Energy = "+str(solution.x[5]))
        self.data = None
        # include both linefit and improved recos for comparison
        frame[self.output] = recoParticle
        frame[self.output + "_nloglike"] = dataclasses.I3Double(results[min_index].fun)
        frame[self.output + "_seed"] = dataclasses.I3String(self.seedtrack[min_index])
        self.PushFrame(frame)
