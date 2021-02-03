'''
Functions commonly used for likelihood analysis
'''
from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import argparse
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi2
import scipy
from tabulate import tabulate

'''

Functions used

'''

def gaussian(x, pos, wid, amp):
    y = amp*np.exp(-4*np.log(2)*((x-pos)/(wid))**2)
    return y

def biGauss(x, pos, wid, r, amp):
    mask = x < pos

    y_all = ([])
    for i in range(0, len(mask)):

        if mask[i] == True:
            m = 1
            nm = 0
        else:
            m = 0
            nm = 1
        if r != 0:
            y1 = gaussian(x[i],pos,r*wid/(r+1),amp)*m
            y2 = gaussian(x[i],pos,wid/(r+1),amp)*nm
            y = y1 + y2
        else:
            y = gaussian(x[i],pos,wid, amp)*nm

        y_all = np.append(y_all, y)
    return y_all

def double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2):
    b1 = biGauss(x, pos1, wid1, r1, amp1)
    b2 = biGauss(x, pos2, wid2, r2, amp2)
    b = np.append(b1, b2)
    return b1+b2

def log_likelihood_biGauss(theta, n, x, debug):
    pos, wid, r, amp = theta
    model = biGauss(x, pos, wid, r, amp)
    L = model - (n*np.log(model))
    if debug == True:
        #print('*****************One Peak*****************')
        print(tabulate([[pos, wid, r, amp, np.sum(L)]], tablefmt=u'fancy_grid',
        headers=("pos", "wid", "r", "amp", "log likelihood")))
    return np.sum(L)

def log_likelihood_doublePeak(theta, n, x, debug):
    pos1, wid1, r1, amp1, pos2, wid2, r2, amp2 = theta
    model = double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2)
    L = model - (n*np.log(model))
    if debug == True:
        #print('*****************Double Peak*****************')
        print(tabulate([[pos1, wid1, r1, amp1, pos2, wid2, r2, amp2, np.sum(L)]], tablefmt=u'fancy_grid',
        headers=("pos1", "wid1", "r1", "amp1", "pos1", "wid1", "r1", "amp1", "log likelihood")))
    return np.sum(L)

def expGauss(x, pos, wid, k, amp):
    aux = (x-pos)/wid
    #val = amp*stats.exponnorm.pdf(aux,k)
    l = 1/(wid*k)
    x_exp = l*(pos - x + (l*wid**2/2))
    x_erf = (pos + l*wid**2 - x)/(np.sqrt(2)*wid)
    val = amp * np.exp(x_exp) * (scipy.special.erfc(x_erf))
    return val

def expDoublePeak(x, pos1, wid1, k1, amp1, pos2, wid2, k2, amp2):
    b1 = expGauss(x, pos1, wid1, k1, amp1)
    b2 = expGauss(x, pos2, wid2, k2, amp2)
    return b1+b2

def log_likelihood_expGauss(theta, n, x, debug):
    pos, wid, k, amp = theta
    model = expGauss(x, pos, wid, k, amp)
    L = model - (n*np.log(model))
    if debug == True:
        #print('*****************One Peak*****************')
        #print("\t".join(["%0.5f" % x for x in [np.sum(L), pos, wid, k, amp]]))
        print(tabulate([[pos, wid, k, amp, np.sum(L)]], tablefmt=u'fancy_grid',
        headers=("pos", "wid", "k", "amp", "log likelihood")))
    return np.sum(L)

def log_likelihood_expDoublePeak(theta, n, x, debug):
    pos1, wid1, k1, amp1, pos2, wid2, k2, amp2 = theta
    model = expDoublePeak(x, pos1, wid1, k1, amp1, pos2, wid2, k2, amp2)
    L = model - (n*np.log(model))
    if debug == True:
        #print('*****************Double Peak*****************')
        #print("\t".join(["%0.5f" % x for x in [np.sum(L), pos1, wid1, k1, amp1, pos2, wid2, k2, amp2]]))
        print(tabulate([[pos1, wid1, k1, amp1, pos2, wid2, k2, amp2, np.sum(L)]], tablefmt=u'fancy_grid',
        headers=("pos1", "wid1", "k1", "amp1", "pos2", "wid2", "k2", "amp2", "log likelihood")))
    return np.sum(L)

def likelihood_ratio_doublePeak(x, n, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2):
    #Likelihood ratio for poisson distributions
    model = double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2)
    val = model - n + (n*np.log(n/model))
    #print('log - ', n/model, 'n - ', n)
    return np.sum(val)

def likelihood_ratio_biGauss(x, n, pos, wid, r, amp):
    #Likelihood ratio for poisson distributions
    model = biGauss(x, pos, wid, r, amp)
    val = model - n + (n*np.log(n/model))
    #print('log - ', n/model, 'n - ', n)
    return np.sum(val)
