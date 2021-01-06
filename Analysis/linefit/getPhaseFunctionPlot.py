#!/usr/bin/env python

import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt

def hgSample(theta, g):
    return (0.5*(1-g**2)/(1+g**2-2*g*cos(np.radians(theta)))**(3/2.0))*np.pi/180

def molecular_scatt(x):
    beta = 0.06225 * (1 + 0.835 * cos(np.radians(x))*cos(np.radians(x)))
    return beta


particle_angle = np.append([0.1,
        0.12589,
        0.15849,
        0.19953,
        0.25119,
        0.31623,
        0.39811,
        0.50119,
        0.63096,
        0.79433,
        1.0000,
        1.2589,
        1.5849,
        1.9953,
        2.5119,
        3.1623,
        3.9811,
        5.0119,
        6.3096,
        7.9433], np.arange(10,185,5))

particle_phase = [1.76661e3,
                          1.29564e3,
                          9.50172e2,
                          6.99092e2,
                          5.13687e2,
                          3.76373e2,
                          2.76318e2,
                          2.18839e2,
                          1.44369e2,
                          1.02241e2,
                          7.16082e1,
                          4.95803e1,
                          3.39511e1,
                          2.28129e1,
                          1.51622e1,
                          1.00154e1,
                          6.57957,
                          4.29530,
                          2.80690,
                          1.81927,
                          1.15257,
                          4.89344e-1,
                          2.44424e-1,
                          1.47151e-1,
                          8.60848e-2,
                          5.93075e-2,
                          4.20985e-2,
                          3.06722e-2,
                          2.27533e-2,
                          1.69904e-2,
                          1.31254e-2,
                          1.04625e-2,
                          8.48826e-3,
                          6.97601e-3,
                          5.84232e-3,
                          4.95306e-3,
                          4.29232e-3,
                          3.78161e-3,
                          3.40405e-3,
                          3.11591e-3,
                          2.91222e-3,
                          2.79696e-3,
                          2.68568e-3,
                          2.57142e-3,
                          2.47603e-3,
                          2.37667e-3,
                          2.32898e-3,
                          2.31308e-3,
                          2.36475e-3,
                          2.50584e-3,
                          2.66183e-3,
                          2.83472e-3,
                          3.03046e-3,
                          3.09206e-3,
                          3.15366e-3]

# Input: scattering_angle
particle_scatt_big_ang = lambda x: np.interp(x,particle_angle,particle_phase)

# If scattering angle is < 0.12589, use:
particle_scatt_small_ang = lambda x: 1.29564e3 * np.power(0.12589/x, 1.346)

# Input: scattering_angle (deg)
def particle_scatt(x):
    phase = []
    for ang in x:
        if ang < 0.12589:
            phase = np.append(phase,particle_scatt_small_ang(ang))
        else:
            phase = np.append(phase,particle_scatt_big_ang(ang))
    return phase

eta = 0.6
g = 0.924*(1-eta)
theta = np.linspace(0.001, 175)
dx = theta[1] - theta[0]

icecube = hgSample(theta,g)*sin(np.radians(theta))
icecube = icecube / sum(icecube*dx)

molecular = molecular_scatt(theta)*sin(np.radians(theta))
molecular = molecular/sum(molecular*dx)
particle = particle_scatt(theta)*sin(np.radians(theta))
particle = particle/sum(particle*dx)
matthew = eta*molecular + (1-eta)*particle

plt.semilogy(theta, matthew, label = 'STRAW data analysis')
plt.semilogy(theta, icecube, label = 'CLSim model')
plt.xlabel('Scattering Angle (degrees)')
plt.ylabel('Normalized Counts')
plt.title('Comparison of CLSim and STRAW Scattering Angle Distributions')
plt.legend()
plt.show()
