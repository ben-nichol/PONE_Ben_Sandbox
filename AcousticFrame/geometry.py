'''
Implements realistic string bending based on simulations by Christian Spannfellner and
interpolations by Felix Henningsen and Hamish Johnson.
https://p-one.atlassian.net/wiki/spaces/SA/pages/238649377/Acoustic+positioning+study+-+2022
'''

from icecube import icetray, dataclasses
import numpy as nn

tilt_coefficients = nn.array(
  [[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
   [ 0.00000000e+00, 1.89675157e-02,-1.88782797e-05, 6.86966811e-09],
   [-4.12837591e-11, 3.08289997e+00,-2.57124952e-03, 7.68681299e-07],
   [ 9.25370891e-14,-2.64904419e+00, 2.89333720e-03,-1.13182809e-06]])

def tilt_delta_x(v_x, z):
'''Calculates the horizontal offset due to currents, based on current speed v_x, and height z'''
    v_x = v_x / icetray.I3Units.m * icetray.I3Units.s
    z   = z / icetray.I3Units.m
    dx  = 0
    for i in range(tilt_coefficients.shape[0]):
        for j in range(tilt_coefficients.shape[1]):
            dx += tilt_coefficients[i,j] * v_x**i * z**j
    return dx * icetray.I3Units.m

def createTiltedGeometry(current_speed=0,
                         current_angle=0,
                         string_positions=[[0,0]]):
'''Takes a current speed, an angle [0, 2pi] of the current, and a list of 2D string positions [[x1,y1],...].
Returns an I3Geometry with tilted strings.'''
    if not 0 <= current_speed <= 0.2 * icetray.I3Units.m / icetray.I3Units.s:
        print('Current speed outside of range for which fit '
              'was performed. Results will not be reliable.'
              'See https://p-one.atlassian.net/wiki/spaces/'
              'SA/pages/238649377/Acoustic+positioning+study+-+2022')

    nominal_module_heights = nn.linspace(50, 1000, 20) * icetray.I3Units.m
    horizontal_offset      = tilt_delta_x(current_speed, nominal_module_heights)
    h  = nn.insert(nominal_module_heights, 0, 0)
    r  = nn.insert(horizontal_offset,      0, 0)
    r[r<0] = 0
    dh = nn.sqrt(nn.diff(h)**2 - nn.diff(r)**2)
    actual_module_heights = nn.cumsum(dh)
    string_tilt           = nn.arctan(nn.diff(r) / nn.diff(h))

    geometry = dataclasses.I3Geometry()
    omgeomap = geometry.omgeo
    for i, (string_x, string_y) in enumerate(string_positions):
        for j, _ in enumerate(actual_module_heights):
            omkey = icetray.OMKey(i + 1, j + 1)
            omdir = dataclasses.I3Direction(
                string_tilt[j], nn.random.uniform(low=0, high=2*nn.pi))
            ompos = dataclasses.I3Position(
                string_x + horizontal_offset[j] * nn.cos(current_angle), 
                string_y + horizontal_offset[j] * nn.sin(current_angle),
                actual_module_heights[j] - 500 * icetray.I3Units.m)
            newomgeo             = dataclasses.I3OMGeo()
            newomgeo.omtype      = dataclasses.I3OMGeo.OMType.IceCube
            newomgeo.orientation = dataclasses.I3Orientation(omdir)
            newomgeo.position    = ompos
            omgeomap[omkey]      = newomgeo
    return geometry

