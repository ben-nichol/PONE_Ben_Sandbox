import numpy as np
from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate
from scipy.optimize import minimize, rosen, rosen_der
import scipy.optimize as op
import pickle

def cpandel(t, d, n = 0.0, sigma = 2.0, lambda_s = 120., rho = 0.004, t0=0.0):

	t -= t0 + d*n

	xi = d/lambda_s                                                             
	eta = rho*sigma - (t/sigma)
	if t<-25.*sigma or t>3500. :
		return 0.0

	if (t>-5.0*sigma and t<30.0*sigma) and xi<5.0 :
		# Define our region dependent approximations of the CPandel function
		_pdf = sp.hyp1f1(0.5*xi,0.5,0.5*eta**2)/sp.gamma(0.5*(xi + 1.))
		_pdf -= np.sqrt(2.)*eta*sp.hyp1f1(0.5*(xi+1.),1.5,0.5*eta**2)/sp.gamma(0.5*xi) 
		_pdf *= (rho**xi)*(sigma**(xi - 1.))*np.exp(-(t**2)/(2.*sigma**2))
		_pdf /= 2.**((1.+xi)/2.)
		return _pdf

	if xi <= 1. and t > 30.*sigma :
		_pdf = np.exp((rho**2)*(sigma**2)/2.)
		_pdf *= (rho**xi)*(t**(xi-1.))*np.exp(-rho*t)
		_pdf /= sp.gamma(xi)
		return _pdf

	if xi>1.0 and t>(rho*(sigma**2.0)) :
		z = max(0.0,-eta/np.sqrt(4*xi - 2.))
		k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
		beta = 0.5*((z/np.sqrt(1. + z**2)) - 1.)
		N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
		N2 = ((beta**2)/(288.))*(6160*beta**4.0 + 18480*beta**3.0 + 19404*beta**2.0 + 8028*beta + 945.)
		phi = 1. - N1/(2.*xi - 1.) + N2/((2.*xi - 1.)**2)
		alpha = -t**2/(2*sigma**2)+0.25*eta**2 - xi*0.5+0.25 + k*(2*xi - 1.) - 0.25*np.log(1 + z**2) - 0.5*xi*np.log(2) + 0.5*(xi-1.)*np.log(2*xi-1.) + xi*np.log(rho) + (xi-1.)*np.log(sigma)
		_pdf = np.exp(alpha)*phi/sp.gamma(xi)
		return _pdf

	if xi>1.0 and t<=(rho*(sigma**2.0)) :
		z = max(0.0,eta/np.sqrt(4*xi-2.))
		k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
		beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
		N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
		N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
		psi = 1. + N1/(2*xi - 1.) + N2/((2*xi - 1.)**2)
		_pdf = (rho**xi)*(sigma**(xi-1.))
		_pdf *= np.exp(0.25*(eta**2.0)-(t**2)/(2*sigma**2))
		_pdf /= np.log(2.0*np.pi)
		U = np.exp(0.5*xi - 0.25)*((2*xi - 1.)**(-0.5*xi))*(2.**(0.5*(xi - 1.)))
		_pdf += U
		_pdf *= np.exp(-k*(2*xi-1.))
		_pdf *= (1. + z**2)**(-0.25)
		_pdf *= psi
		return _pdf

	if xi<=1. and t<=(rho*(sigma**2.0)) :
		_pdf = (rho*sigma)**xi
		_pdf *= eta**(-xi)
		_pdf *= np.exp(-t**2.0/(2.0*sigma**2.0))
		_pdf /= np.sqrt(2.*np.pi*sigma**2.0)
		return _pdf 

	return 0.0

def readTables() :

    pdf = []
    peaktime = []

    for i in range(400) :
        pdf.append([])
        for j in range(10010) :
            pdf[-1].append(0.0)

    input_dict  = pickle.load(open("ChargeTimePdf_April19th.pkl", 'rb'))

    dist  = input_dict["dist"]
    char  = input_dict["char"]
    time  = input_dict["time"]

    for i in range(len(char)):
        pdf[dist[i]][time[i]] = char[i]

    for i in range(len(pdf)) :
        peaktime.append(0)
        for j in range(len(pdf[i])):
            if pdf[i][i] > pdf[i][peaktime[-1]]:
                peaktime[-1] = j

    return pdf, peaktime

def LikelihoodFunctor(_pdf,_peaktime):

    pdf = _pdf
    peaktime = _peaktime

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(n=0., sigma = 2.0, lambda_s = 120., rho = 0.004, t0=0.0, atten=50.,_amp=1.0):
        chi2 = 0.0;
        #for d in range(10,11) :
        for d in range(10,100) :
                for t in range(-5,100) :
                        mc = pdf[d][t+10]
                        amp = _amp*np.exp(-d/atten)/(d*d)
                        amp *= pdf[d][peaktime[d]]/cpandel(float(peaktime[d])-9.5, float(d)+0.5, n, sigma, lambda_s, rho,t0)
                        pandel = amp*cpandel(float(t)+0.5, float(d)+0.5, n, sigma, lambda_s, rho,t0)
                        if mc > 1.0e-8 :
                                chi2 += ((mc-pandel)**2.0)/mc

        return chi2

    return likelihoodFunction

if __name__ == "__main__":

    pdf, peaktime = readTables()

    print(pdf[10])

    qFunctor = LikelihoodFunctor(pdf, peaktime)

    def func(x):
        n,sigma,lambda_s,rho,t0,atten,amp = x
        return qFunctor(n,sigma,lambda_s,rho,t0,atten,amp)

    solution = op.minimize(fun=func,
                               x0=np.array([0.,2.0,120.,0.004,0.0,50.,1.0]),
                               method='Nelder-Mead')

    print(solution)
