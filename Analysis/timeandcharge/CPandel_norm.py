import numpy as np
from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate

def cpandel(t, d, n = 0.0, sigma = 1.1339139328144132, lambda_s = 317.50178764954626, rho = 0.04079084329979382, t0=0.0):

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

pdf = [[]]
time_lim = [0.0,0.0]
dist_lim = [0.0,0.0]
peaktime = [0.0]

def readTables() :
	global pdf
	global time_lim
	global dist_lim
	global peaktime

	infile = open("fittertables.dat","r")
	lines = infile.readlines()
	linecount = 0
	xcount = 0
	ny = 0
	nx = 0
	minx = 0.0
	maxx = 0.0
	miny = 0.0
	maxy = 0.0
	maxvalue = 0.0;

	for line in lines :
		splitline = line.split(",",100)
		if linecount == 0 :
			nx = int(splitline[0].replace("\n",""))
			ny = int(splitline[1].replace("\n",""))
			minx = float(splitline[2].replace("\n",""))
			maxx = float(splitline[3].replace("\n",""))
			miny = float(splitline[4].replace("\n",""))
			maxy = float(splitline[5].replace("\n",""))
			linecount += 1
		else :
			if xcount == ny :
				pdf.append([])
				maxvalue = 0.0
				peaktime.append(0.0)
				xcount = 0
			for value in splitline :
				pdf[-1].append(float(value.replace("\n","")))
				if pdf[-1][-1] > maxvalue :
					maxvalue = pdf[-1][-1]
					peaktime[-1] = len(pdf[-1])-1
				xcount += 1
			linecount += 1
	time_lim = [miny,maxy]
	dist_lim = [minx,maxx]
	print(len(pdf))

def ComputeChiSqr(n=0., sigma = 2.0, lambda_s = 120., rho = 0.004, t0=0.0) :
	global pdf
	global peaktime

	chi2 = 0.0;
	#for d in range(10,11) :
	for d in range(10,100) :
		for t in range(-5,100) :
			mc = pdf[d][t+10]
			amp = pdf[d][peaktime[d]]/cpandel(float(peaktime[d])-9.5, float(d)+0.5, n, sigma, lambda_s, rho,t0)
			pandel = amp*cpandel(float(t)+0.5, float(d)+0.5, n, sigma, lambda_s, rho,t0)
			if mc > 1.0e-8 :
				chi2 += ((mc-pandel)**2.0)/mc

	return chi2

if __name__ == "__main__":

	d_10=[]
	d_20=[]
	d_50=[]
	d_100=[]
	t = []

	for i in range(10010) :
		d_10.append(cpandel(-10.+float(i),10.0))
		d_20.append(cpandel(-10.+float(i),20.0))
		d_50.append(cpandel(-10.+float(i),50.0))
		d_100.append(cpandel(-10.+float(i),100.0))

	integral_10 = sum(d_10)
	integral_20 = sum(d_20)
	integral_50 = sum(d_50)
	integral_100 = sum(d_100)

	print("integral_10 = "+str(integral_10)+" integral_20 = "+str(integral_20)+" integral_50 = "+str(integral_50)+ " integral_100 = "+str(integral_100))  


