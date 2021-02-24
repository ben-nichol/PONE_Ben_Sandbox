'''
A single place to keep all PDFs that have been coded in the process of attempting to create a muon track reconstruction for P-ONE. 

'''

import numpy as np
import math

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate
import sys
import inspect

# Some quantities that are environment dependent
c = 2.99792458e8                                # speed of light 
n = 1.34                                        # 1.33 is the refractive index of water at 20 degrees C
c_m = c/n                                       # light in water
theta_c = np.arccos(1./n)                       # Cherenkov angle in water
#lambda_s = 120 * I3Units.m                     # scattering length of light for violet light
#lambda_a = 15 * I3Units.m                      # absorption length of light for violet light
#tau = 557E-9 * I3Units.second                  # time parameter that has to be fit using simulations or data      

# Pandel function
def pandel(t, d, lambda_a = 15., lambda_s = 120., tau = 557E-9):
    N = np.exp(-d/lambda_a)*np.power((1. + (tau*c_m)/lambda_a),-d/lambda_s)
    exp = np.exp(-t*((1./tau)+(c_m/lambda_a))-d/lambda_a)
    frac = (np.power(tau,-d/lambda_s)*np.power(t,(d/lambda_s) - 1))/(sp.gamma(d/lambda_s))
    return (1./N)*frac*exp

# CPandel function. This only works on a constrained subset of the t-d plane, namely for direct hits (time is in ns) #lambda_s = 33.3 (default)
#def cpandel(t, d, sigma = 10., lambda_s = 120., rho = 0.004):
#    xi = d/(lambda_s*np.sin(theta_c))
#    eta = rho*sigma - (t/sigma)
#    frac = (np.power(rho,xi)*np.power(sigma,xi-1.)*np.exp(-(t**2)/(2.*sigma**2)))/(np.power(2.,(1. + xi)/2.))
#    frac_1 = sp.hyp1f1(0.5*xi,0.5,0.5*eta**2)/sp.gamma(0.5*(xi + 1.))
#    frac_2 = sp.hyp1f1(0.5*(xi+1),1.5,0.5*eta**2)/sp.gamma(0.5*xi)
#    return frac*(frac_1 - np.sqrt(2.)*eta*frac_2)

def cpandel(t, d, sigma = 1.1339139328144132, lambda_s = 317.50178764954626, rho = 0.04079084329979382):

    pdf = []
    for i in range(len(t)) :
        xi = d[i]/lambda_s                                                             
        eta = rho*sigma - (t[i]/sigma)
        if t[i]<-25.*sigma or t[i]>3500. :
            pdf.append(0.0)

        elif (t[i]>-5.0*sigma and t[i]<30.0*sigma) and xi<5.0 :
            # Define our region dependent approximations of the CPandel function
            _pdf = sp.hyp1f1(0.5*xi,0.5,0.5*eta**2)/sp.gamma(0.5*(xi + 1.))
            _pdf -= np.sqrt(2.)*eta*sp.hyp1f1(0.5*(xi+1.),1.5,0.5*eta**2)/sp.gamma(0.5*xi) 
            _pdf *= (rho**xi)*(sigma**(xi - 1.))*np.exp(-(t[i]**2)/(2.*sigma**2))
            _pdf /= 2.**((1.+xi)/2.)
            pdf.append(_pdf)

        elif xi <= 1. and t[i] > 30.*sigma :
            _pdf = np.exp((rho**2)*(sigma**2)/2.)
            _pdf *= (rho**xi)*(t[i]**(xi-1.))*np.exp(-rho*t[i])
            _pdf /= sp.gamma(xi)
            pdf.append(_pdf)

        elif xi>1.0 and t[i]>(rho*(sigma**2.0)) :
            z = max(0.0,-eta/np.sqrt(4*xi - 2.))
            k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
            beta = 0.5*((z/np.sqrt(1. + z**2)) - 1.)
            N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
            N2 = ((beta**2)/(288.))*(6160*beta**4.0 + 18480*beta**3.0 + 19404*beta**2.0 + 8028*beta + 945.)
            phi = 1. - N1/(2.*xi - 1.) + N2/((2.*xi - 1.)**2)
            alpha = -t[i]**2/(2*sigma**2) + 0.25*eta**2 - xi*0.5 + 0.25 + k*(2*xi - 1.) - 0.25*np.log(1 + z**2) - 0.5*xi*np.log(2) + 0.5*(xi-1.)*np.log(2*xi-1.) + xi*np.log(rho) + (xi-1.)*np.log(sigma)
            _pdf = np.exp(alpha)*phi/sp.gamma(xi)
            pdf.append(_pdf)

        elif xi>1.0 and t[i]<=(rho*(sigma**2.0)) :
            z = max(0.0,eta/np.sqrt(4*xi-2.))
            k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
            beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
            N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
            N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
            psi = 1. + N1/(2*xi - 1.) + N2/((2*xi - 1.)**2)
            _pdf = (rho**xi)*(sigma**(xi-1.))*np.exp(0.25*(eta**2.0)-(t[i]**2)/(2*sigma**2))
            _pdf /= np.log(2.0*np.pi)
            U = np.exp(0.5*xi - 0.25)*((2*xi - 1.)**(-0.5*xi))*(2.**(0.5*(xi - 1.)))
            _pdf += U
            _pdf *= np.exp(-k*(2*xi-1.))
            _pdf *= (1. + z**2)**(-0.25)
            _pdf *= psi
            pdf.append(_pdf)

        elif xi<=1. and t[i]<=(rho*(sigma**2.0)) :
            _pdf = (rho*sigma)**xi
            _pdf *= eta**(-xi)
            _pdf *= np.exp(-t[i]**2.0/(2.0*sigma**2.0))
            _pdf /= np.sqrt(2.*np.pi*sigma**2.0)
            pdf.append(_pdf) 

    return pdf

# -log(CPandel) so that the result can be a sum instead of a product. Multiple cases are used as approximations are needed for different domains of t-d. #
# This region restriction has made the python implementation more complicated in the goals of keeping it pythonic (and hence fast-ish). In particular,   #
# this function builds other functions that are the possible cases, then finds which elements of the numpy arrays satisfy those particular conditionals  #
# and runs those particular arrays through it's respective algorithm. At the end we rebuild our output array before returning.
# time is in nanoseconds

# NOTE: Could Potentially be bugged (lambda_s = 33.3)
def log_cpandel(t, d, sigma = 10, lambda_s = 120., rho = 0.004):
    xi = d/lambda_s
    eta = rho*sigma - (t/sigma)
    darkprob = 1./10000.   

    # Define our region dependent approximations of the CPandel function
    def region1(time, dist, eta_in, xi_in):
        first = xi_in*np.log(rho) + (xi_in - 1.)*np.log(sigma) - (time**2)/(2.*sigma**2) - ((1. + xi_in)/2.)*np.log(2)
        frac_1 = sp.hyp1f1(0.5*xi_in,0.5,0.5*eta_in**2)/sp.gamma(0.5*(xi_in + 1.))
        frac_2 = sp.hyp1f1(0.5*(xi_in+1.),1.5,0.5*eta_in**2)/sp.gamma(0.5*xi_in)
        second = np.log(frac_1 - np.sqrt(2.)*eta_in*frac_2)
        return -(first + second) + darkprob

    def region2(time, dist, eta_in, xi_in):
        first = (rho**2)*(sigma**2)/2.
        second = xi_in*np.log(rho) + (xi_in - 1.)*np.log(time) - np.log(xi_in) - rho*time - sp.gamma(xi_in)
        return -(first + second) + darkprob   

    def region3(time, dist, eta_in, xi_in):
        #print("r3 Avg time = " + str(np.average(time)) + " and average d = " + str(np.average(dist)))
        #print("avg xi = " + str(np.average(xi_in)))
        z = -eta_in/np.sqrt(4*xi_in - 2.)
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
        phi = 1. - N1/(2.*xi_in - 1.) + N2/((2.*xi_in - 1.)**2)
        alpha = -time**2/(2*sigma**2) + 0.25*eta_in**2 - xi_in*0.5 + 0.25 + k*(2*xi_in - 1.) - 0.25*np.log(1 + z**2) - 0.5*xi_in*np.log(2) + 0.5*(xi_in - 1.)*np.log(2*xi_in - 1.) + xi_in*np.log(rho) + (xi_in - 1.)*np.log(sigma)
        return -(alpha - np.log(sp.gamma(xi_in)) + np.log(phi)) + darkprob   

    def region4(time, dist, eta_in, xi_in):
        #print("r4 Avg time = " + str(np.average(time)) + " and average d = " + str(np.average(dist)))
        z = eta_in/np.sqrt(4*xi_in - 2.)
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
        psi = 1. + N1/(2*xi_in - 1.) + N2/((2*xi_in - 1.)**2)
        first = xi_in*np.log(rho) + (xi_in - 1.)*np.log(sigma) - (time**2)/(2*sigma**2) + 0.25*(eta_in**2) - 0.25*np.log(2*np.pi)
        U = 0.5*xi_in - 0.25 - 0.5*xi_in*np.log(2*xi_in - 1.) + 0.5*(xi_in - 1.)*np.log(2)
        second = -k*(2*xi_in - 1.) - 0.25*np.log(1. + z**2) + np.log(psi)
        return -(first + U + second) + darkprob  

    def region5(time, dist, eta_in, xi_in):
        return -(xi_in*np.log(rho*sigma) - 0.5*np.log(2*np.pi*sigma**2) -
            xi_in*np.log(eta_in) - (time**2)/(2*sigma**2)) + darkprob   
        
    # Now we find all the indices corresponding with their respective regions

    # The case where we can use the exact form since t and d are small enough
    r1 = np.where(np.logical_and(np.less_equal(xi, 5.), np.logical_and(np.less_equal(rho*sigma - (300)/sigma, eta), np.less(eta, rho*sigma + 5.)))) 
    # "Small" distance, but large and positive t
    r2 = np.where(np.logical_and(np.less_equal(xi, 1.), np.less(eta, rho*sigma - (300)/sigma)))
    # Large distance and large + positive t
    r3 = np.where(np.logical_or(np.logical_and(np.greater(xi, 5.), np.less(eta, 0.)),
                                np.logical_and(np.greater(xi, 1.), np.less(eta, rho*sigma - (300)/sigma))))
    # Large distance and "small" t
    r4 = np.where(np.logical_or(np.logical_and(np.greater(xi, 5.), np.greater_equal(eta, 0.)),
                                np.logical_and(np.greater(xi, 1.), np.greater_equal(eta, rho*sigma + 5.))))
    # small distance and negative time
    r5 = np.where(np.logical_and(np.less_equal(xi, 1.), np.greater_equal(eta, rho*sigma + 5.)))

    # check if everything is satisfied
    test = len(r1[0]) + len(r2[0]) + len(r3[0]) + len(r4[0]) + len(r5[0])
    if test != len(t):
        print("Invalid domain recieved. Oh no. Values are (t,d) = " + str((t,d)))
        print("Array lengths summed to " + str(test))
        print("Total length should be " + str(len(t)))
        if math.isnan(t[0]):
            print("Nan found. Returning Nan")
        else:
            sys.exit(0)

    # Run the stuff
    
    result1 = region1(t[r1], d[r1], eta[r1], xi[r1])
    result2 = region2(t[r2], d[r2], eta[r2], xi[r2])
    result3 = region3(t[r3], d[r3], eta[r3], xi[r3])
    result4 = region4(t[r4], d[r4], eta[r4], xi[r4]) 
    result5 = region5(t[r5], d[r5], eta[r5], xi[r5])

    # Debugging
    
    #print("r1 = " + str(t[r1]))
    #print("r2 = " + str(t[r2]))
    #print("r3 = " + str(t[r3]))
    #print("r4 = " + str(t[r4]))
    #print("r5 = " + str(t[r5]))
    
#    print("r3 sum = " + str(np.sum(result3))) 
#    print("r4 sum = " + str(np.sum(result4))) 

    # Replace result in correct locations and return 
    result = np.zeros(len(t))
    
    result[r1] = result1
    result[r2] = result2
    result[r3] = result3  
    result[r4] = result4  
    result[r5] = result5

    if math.isnan(t[0]):
        result = t

    return result

# Reparamaterized Pandel used in derivation of CPandel
def pandel2(t, d, lambda_s = 120., rho = 0.004):
    xi = d/(lambda_s*np.sin(theta_c))    
    num = np.power(rho, xi)*np.power(t,xi - 1)
    denom = sp.gamma(xi)
    exp = np.exp(-rho*t)
    return (num/denom)*exp

# A test PDF that is not distance dependent
# Modified CPandel function that is fixed in distance and extended for all time (time is in ns) #lambda_s = 33.3 (default)
def mod_logcpandel(t, d, sigma = 10., lambda_s = 120., rho = 0.004):
    d = 40.*np.ones(len(t))           # Fix the length for our purposes
    xi = d/lambda_s                   # *np.sin(theta_c) but this is accounted for in initial distance computations
    eta = rho*sigma - (t/sigma)
    
    def region1(time, dist, xi_in, eta_in):
        return -np.log(cpandel(time, dist, sigma, lambda_s, rho))

    def region2(time, dist, xi_in, eta_in):
        first = np.exp((rho**2)*(sigma**2)/2.)
        second = pandel2(time, dist, lambda_s, rho)
        return -np.log(first*second)

    def region5(time, dist, xi_in, eta_in):
        num = np.power(rho*sigma,xi_in)
        denom = np.sqrt(2*np.pi*sigma**2)
        prod1 = np.power(eta_in,-xi_in)
        exp = np.exp(-(time**2)/(2.*sigma**2))
        return -np.log((num/denom)*prod1*exp)

    # Now we find when to switch to each approximation

    # exact region -50 < t < 300
    r1 = np.where(np.logical_and(np.less_equal(t, 300.), np.greater_equal(t, -5.*sigma)))
    # large t > 300
    r2 = np.where(np.greater(t, 300.))
    # small t < -50
    r5 = np.where(np.logical_and(np.less(t, -5.*sigma), np.greater_equal(t, -360)))
    r5_pre = np.where(np.less(t, -360))

    # check if everything is satisfied
    test = len(r1[0]) + len(r2[0]) + len(r5[0]) + len(r5_pre[0])
    if test != len(t):
        print(r1[0])
        print(r2[0])
        print(r5[0])
        print(r5_pre[0])
        print("Invalid domain recieved. Oh no. Values are (t,d) = " + str((t,d)))
        print("Array lengths summed to " + str(test))
        print("Total length should be " + str(len(t)))
        if math.isnan(t[0]):
            print("Nan found. Returning Nan")
        else:
            sys.exit(0)

    # compute the actual value
    result1 = region1(t[r1], d[r1], xi[r1], eta[r1])
    result2 = region2(t[r2], d[r2], xi[r2], eta[r2])
    result5 = region5(t[r5], d[r5], xi[r5], eta[r5])
    
    # special regions
    t_pre = -np.ones(len(r5_pre[0]))*360.
    eta_pre = rho*sigma - (t_pre/sigma)
    result5_pre = region5(t_pre, d[r5_pre], xi[r5_pre], eta_pre)


        # Build final result    
    result = np.zeros(len(t))

    result[r1] = result1
    result[r2] = result2
    result[r5] = result5
    result[r5_pre] = result5_pre

    if math.isnan(t[0]):
        result = t

    return result

# Using an exponentially modified gaussian to build a sharper distribution 

#def xnorm_pdf(t, d, sigma = 10., lambda_s = 120., rho = 0.004):




# Use interpolate and produce a pdf based off of a particular data set from STRAW. 
# Ideally we would normalize via integration, and find the 'zero' residual time by integrating and comparing with cpandel. 
# For now, however, it is simple and easier to normalize by matching the peak value with the cpandel peak. Then the zero is matched the same way. 

# First we import the data needed
straw = np.load("/home/users/ghuman/simAnalysis/output/npy/20200630_201044_UTC_SDOM1_FLASHRUN_EDGECOMBE_FLASH_P1_both_1000Hz_blue_8V_RUN1_120s_2020-06-30_2011_20182201054_ch1_cpan-normalized.npy")

# ../straw-b/pdf_plot.py tells us that the 'zero' bin is 434543 when compared with cpandel(t, 40). So, we can build our domain using this 
# Note: Peak is at bin 303 in the following regime
x = np.array(range(-300, 3000))
y = straw[(x[0] + 434543):(x[-1] + 434544)]

# integrate using simpsons rule to normalize
#print(integrate.simps(y)) # = 0.281694727151
y = np.divide(y, 0.281694727151)

# Interpolated function using linear fit
f = inter.interp1d(x, y)

# Interpolated function using a cubic method
f2 = inter.interp1d(x, y, kind='cubic')

# splinefit
f3 = inter.splrep(x, y)

# Removing zero values and setting to the next lowest value that occurs, since I will want to eventually compute the logarithm

# Find what value is the second lowest
above_zero = np.where(np.greater(y, 0.0))
minimum = min(y[above_zero])

# Find zero values and set to new minimum
below_min = np.where(np.less(y, minimum))
y[below_min] = minimum*np.ones(len(below_min[0]))

# Build a filtered distribution

# use -log of y
ylog = -np.log(y)

# Ken's Work: ---

# Gaussian Filter:
def gaussWindow(xData,yData,stdDev):
    gaussSmoothedVals = np.zeros(yData.size)#set up the output array
    ct=0
    for xPosition in xData:                           #loop over all points
        #print(xPosition)
        kernel = np.exp(-(xData - xPosition)**2/(2*stdDev**2)) #set up the gaussian
        kernel = kernel / sum(kernel)                 #scale the gaussian to preserve the total
        gaussSmoothedVals[ct] = sum(yData * kernel)  #do the multiplication
        ct+=1
    return gaussSmoothedVals                          #return the points

highEdge = 33 #the upper limit of the part which is left alone (the peak)
lowEdge = -12 #the lower limit of the same as above
stdDev = 50 #the standard deviation of the gaussian used to do the smoothing

smoothPointsX = x[x>highEdge] #get the points after the peak to smooth
smoothPointsY = ylog[x>highEdge]

prePointsX = x[x<=lowEdge] #get the points before the peak to smooth
prePointsY = ylog[x<=lowEdge]

okPointsX = x[(x>lowEdge) & (x<=highEdge)] #leave the peak alone
okPointsY = ylog[(x>lowEdge) & (x<=highEdge)] 

gaussSmoothedBegin = gaussWindow(prePointsX,prePointsY,stdDev) #same on the low side
gaussSmoothedMid = gaussWindow(okPointsX,okPointsY,stdDev) # and the peak
gaussSmoothedEnd = gaussWindow(smoothPointsX,smoothPointsY,stdDev) #do the smoothing on the high side

gauss_combined = np.concatenate((gaussSmoothedBegin,gaussSmoothedMid,gaussSmoothedEnd),axis=None) #put the points back together
# ---

savgolBegin = savgol_filter(prePointsY, 101, 4)
savgolMid = savgol_filter(okPointsY, 11, 4)
savgolEnd = savgol_filter(smoothPointsY, 451, 4)

sav_combined = np.concatenate((savgolBegin,savgolMid,savgolEnd),axis=None) #put the points back together

# yfilt = savgol_filter(-np.log(y), 451, 4) # savgol_filter method (not currently used)

f4 = inter.interp1d(x, gauss_combined, kind='cubic')

# Save as numpy array
#xylog = np.zeros((2, len(x)))
#xylog[0,:] = x
#xylog[1,:] = -np.log(y)
#np.save("/home/users/ghuman/simAnalysis/output/npy/log_data_array.npy", xylog)

def straw_pdf(t,d,fit='filter'):
    # If the time is greater than t_max seconds or less than t_min, we can assume it is the minimum value above zero
    time = t.copy()
    t_max = 150.
    t_min = -12.
    too_high = np.where(np.greater_equal(time, t_max))
    too_low = np.where(np.less_equal(time, t_min))
    time[too_high] = t_max*np.ones(len(too_high[0]))
    time[too_low] = t_min*np.ones(len(too_low[0]))
    if fit=='spline':
        return -np.log(inter.splev(time, f3))
    elif fit=='filter':
        return f4(time)
    else:
        return -np.log(f2(time))

'''
rv = stats.exponnorm(K=1, loc=3)
xnew = np.linspace(-50, 100, 10000)

# plot
plt.plot(x, y, linestyle='', marker='.', color='b', label='straw data')
plt.plot(xnew, rv.pdf(xnew), linestyle='-', color='k', label='exponnorm')
plt.title('Exponnorm')
plt.xlabel(r'time(ns)')
plt.ylabel(r'Effective Likelihood')
plt.grid()
plt.legend()
plt.xlim(-25, 50)
plt.savefig('/home/users/ghuman/simAnalysis/output/plots/straw-b/straw_exponnorm.png', dpi=300)
#plt.plot(xnew, f(xnew), linestyle='-', color='r', label='linear interpolated')
#plt.plot(xnew, f2(xnew), linestyle='--', color='g', label='cubic interpolation')
#plt.plot(xnew, ynew, linestyle=':', color='k', label='spline fit')
#plt.plot(x, sav_combined, linestyle='-', color='c', label='savgol filtered')
#plt.plot(x, gauss_combined, linestyle='-', color='r', label='gauss filtered')
#plt.xlim(-50, 150)
#plt.legend()
#plt.savefig('/home/users/ghuman/simAnalysis/output/plots/straw-b/pdf_spline_straw.png', dpi=400)
'''

'''
# plot
plt.plot(x, -np.log(y), linestyle='', marker='.', color='b', label='straw data')
plt.title('STRAW pdf')
plt.xlabel(r'time(ns)')
plt.ylabel(r'Effective Likelihood')
plt.grid()
#plt.savefig('/home/users/ghuman/simAnalysis/output/plots/straw-b/plot_straw.png', dpi=300)
#plt.plot(xnew, f(xnew), linestyle='-', color='r', label='linear interpolated')
#plt.plot(xnew, f2(xnew), linestyle='--', color='g', label='cubic interpolation')
#plt.plot(xnew, ynew, linestyle=':', color='k', label='spline fit')
plt.plot(x, sav_combined, linestyle='-', color='c', label='savgol filtered')
plt.plot(x, gauss_combined, linestyle='-', color='r', label='gauss filtered')
plt.xlim(-50, 150)
plt.legend()
#plt.savefig('/home/users/ghuman/simAnalysis/output/plots/straw-b/pdf_spline_straw.png', dpi=400)
plt.savefig('/home/users/ghuman/simAnalysis/output/plots/straw-b/pdf_interpolation_straw_gauss_zoomed.png', dpi=400)
'''

