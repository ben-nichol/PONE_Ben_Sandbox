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

'''

Functions used

'''

def log_likelihood_biGauss_functor(time,charge):

    timearray = time
    chargearray = charge

    def gaussian(x,mean,sigma):
        return (1./(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mean)**2.0)/(2.*sigma**2))

    def biGauss(x, mean, sigma1,sigma2):
        if x > mean :
            return gaussian(x[i],mean,sigma1)
        else :
            return (sigma1/sigma2)*gaussian(x[i],mean,sigma2)

    def log_likelihood_biGauss(mean,sigma1,sigma2):
        sumloglike = 0.0
        for i in range(len(timearray)) :
            model = biGauss(timearray[i],mean,sigma1,sigma2)
            sumloglike += chargearray[i]*np.log(model)
        return sumloglike

    return log_likelihood_biGauss

def log_likelihood_doublePeak_functor(time,charge):

    timearray = time
    chargearray = charge

    def gaussian(x,mean,sigma):
        return (1./(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mean)**2.0)/(2.*sigma**2))

    def biGauss(x, mean, sigma1,sigma2):
        if x > mean :
            return gaussian(x[i],mean,sigma1)
        else :
            return (sigma1/sigma2)*gaussian(x[i],mean,sigma2)

    def double_peak(mean1,sigma1,mean2,sigma2,r):
        b1 = biGauss(x,mean1,sigma1,sigma2)
        b2 = biGauss(x,mean2,sigma3,sigma4)
        b = np.append(b1, b2)
        return b1+b2

    def log_likelihood_doublePeak(theta, n, x, debug):
        pos1, wid1, r1, amp1, pos2, wid2, r2, amp2 = theta
        model = double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2)
        L = model - (n*np.log(model))
        return np.sum(L)

    return log_likelihood_doublePeak

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
    return np.sum(L)

def log_likelihood_expDoublePeak(theta, n, x, debug):
    pos1, wid1, k1, amp1, pos2, wid2, k2, amp2 = theta
    model = expDoublePeak(x, pos1, wid1, k1, amp1, pos2, wid2, k2, amp2)
    L = model - (n*np.log(model))
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
