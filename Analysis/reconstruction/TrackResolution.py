from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame 
from icecube.dataclasses import I3Particle 
import numpy as np                 
import sys, os
import argparse
import math as m
import ROOT
from scipy import special as sp
from scipy import interpolate as inter
from scipy.signal import savgol_filter
from scipy import stats
from scipy import integrate

PMTacceptance = list()
domaccFile = open("/home/users/tmcelroy/pone_offline/data/config_13.txt","r")
lines = domaccFile.readlines()

zenithcount = 0
for line in lines :
    splitline = line.split(" ",1000)
    if zenithcount % 179 == 0 :
        zenithcount = 0
        PMTacceptance.append([])
    PMTacceptance[-1].append([])
    for value in splitline :
        PMTacceptance[-1][-1].append(float(value))
    zenithcount += 1

def GetPMTAcceptance(photonDir,pmtid):
    global PMTacceptance
 
    if pmtid > len(PMTacceptance)-1 :
         print("Bad PMTID")
         print(pmtid)

    thetaBin = max(0,min(178,int(180.0*photonDir.zenith/np.pi)))
    phiBin = max(0,min(358,int(180.0*photonDir.azimuth/np.pi)))
    if thetaBin > len(PMTacceptance[pmtid]) -1 :
        print("Bad thetaBin")
        print(thetaBin)
    if phiBin > len(PMTacceptance[pmtid][thetaBin]) -1 :
        print("Bad thetaBin")
        print(phiBin)
    return PMTacceptance[pmtid][thetaBin][phiBin]

def cpandel(t, d, sigma = 3.1339139328144132, lambda_s = 317.50178764954626, rho = 0.04079084329979382):

    pdf = 0.0
    xi = d/lambda_s
    scale = 1.0
    if t< -3.0*sigma:
        scale = 1./(1.0+np.abs(t+1))
    if t>3500.0 :
        scale = 1./(1.0+np.abs(t-3500.))
    t = min(3500.,max(-3.0*sigma,t))
    eta = rho*sigma - (t/sigma)

    if (t>-5.0*sigma and t<30.0*sigma) and xi<5.0 :
        # Define our region dependent approximations of the CPandel function
        _pdf = sp.hyp1f1(0.5*xi,0.5,0.5*eta**2)/sp.gamma(0.5*(xi + 1.))
        _pdf -= np.sqrt(2.)*eta*sp.hyp1f1(0.5*(xi+1.),1.5,0.5*eta**2)/sp.gamma(0.5*xi)
        _pdf *= (rho**xi)*(sigma**(xi - 1.))*np.exp(-(t**2)/(2.*sigma**2))
        _pdf /= 2.**((1.+xi)/2.)
        return _pdf*scale

    elif xi <= 1. and t > 30.*sigma :
        _pdf = np.exp((rho**2)*(sigma**2)/2.)
        _pdf *= (rho**xi)*(t**(xi-1.))*np.exp(-rho*t)
        _pdf /= sp.gamma(xi)
        return _pdf*scale

    elif xi>1.0 and t>(rho*(sigma**2.0)) :
        z = max(0.0,-eta/np.sqrt(4*xi - 2.))
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/np.sqrt(1. + z**2)) - 1.)
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4.0 + 18480*beta**3.0 + 19404*beta**2.0 + 8028*beta + 945.)
        phi = 1. - N1/(2.*xi - 1.) + N2/((2.*xi - 1.)**2)
        alpha = -t**2/(2*sigma**2) + 0.25*eta**2 - xi*0.5 + 0.25 + k*(2*xi - 1.) - 0.25*np.log(1 + z**2) - 0.5*xi*np.log(2) + 0.5*(xi-1.)*np.log(2*xi-1.) + xi*np.log(rho) + (xi-1.)*np.log(sigma)
        _pdf = np.exp(alpha)*phi/sp.gamma(xi)
        return _pdf*scale

    elif xi>1.0 and t<=(rho*(sigma**2.0)) :
        z = max(0.0,eta/np.sqrt(4*xi-2.))
        k = 0.5*(z*np.sqrt(1. + z**2) + np.log(z + np.sqrt(1. + z**2)))
        beta = 0.5*((z/(np.sqrt(1. + z**2)) - 1.))
        N1 = (beta/12.)*(20*beta**2 + 30*beta + 9.)
        N2 = ((beta**2)/(288.))*(6160*beta**4 + 18480*beta**3 + 19404*beta**2 + 8028*beta + 945.)
        psi = 1. + N1/(2*xi - 1.) + N2/((2*xi - 1.)**2)
        _pdf = (rho**xi)*(sigma**(xi-1.))*np.exp(0.25*(eta**2.0)-(t**2)/(2*sigma**2))
        _pdf /= np.log(2.0*np.pi)
        U = np.exp(0.5*xi - 0.25)*((2*xi - 1.)**(-0.5*xi))*(2.**(0.5*(xi - 1.)))
        _pdf += U
        _pdf *= np.exp(-k*(2*xi-1.))
        _pdf *= (1. + z**2)**(-0.25)
        _pdf *= psi
        return _pdf*scale

    elif xi<=1. and t<=(rho*(sigma**2.0)) :
        _pdf = (rho*sigma)**xi
        _pdf *= eta**(-xi)
        _pdf *= np.exp(-t**2.0/(2.0*sigma**2.0))
        _pdf /= np.sqrt(2.*np.pi*sigma**2.0)
        return _pdf*scale

    return 0.0

def GetGeoTime(position,vert,direction) :
	c = 0.299792458                                 # speed of light 
	n = 1.34
	ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
	c_n = c/ngroup                                     # light in water
	theta_c = np.arccos(1./n)
	x = position.x - vert.x
	y = position.y - vert.y
	z = position.z - vert.z
	dotprod = x*direction.x + y*direction.y + z*direction.z
        dc = 0.25
        if (x*x + y*y + z*z-dotprod*dotprod) > 0.5 :
	   dc = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
	d = dc/np.sin(theta_c)
        perp = dc/np.tan(theta_c)-dotprod
        pdir = dataclasses.I3Direction(direction.x*perp+position.x,direction.y*perp+position.y,direction.z*perp+position.z)
        t = d/c_n - perp/c
        return d,dc,t,pdir

def timeresidualvariance(vertex,direction,pulse_series,geo_doms):
    dark = 1.e-8
    tau = 18.949132224466762
    sigma = 1.0
    d_o = 50.0
    a_o = 100.00
    d_1 = 0.25

    sum_nloglike = 0.0
    npulses = 0.0
    ndoms = 0.0
    sumcharge = 0.0
    sum_timeres = 0.0
    sum_chargefrac = 0.0
    chargefracs = []
    for dom in pulse_series.keys() :
        domkey =  OMKey(dom.string, dom.om, 0)
        d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
        totalcharge = 0.0
        for pulse in pulse_series[dom] :
            charge = 1.0
            if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                 charge = pulse.charge
            totalcharge += charge
            npulses += 1.0
        A = totalcharge*a_o/np.sqrt(a_o**2.0+totalcharge**2.0)
        D = np.sqrt(d_1**2.0+dc**2.0)
        sum_chargefrac += A*D
        chargefracs.append(A*D)
        sumcharge += A
        ndoms += 1.0

    meanchargefrac = sum_chargefrac/ndoms
    sigmachargefrac = 0.0
    for i in range(len(chargefracs)) :
        sigmachargefrac += (chargefracs[i]-meanchargefrac)**2.0

    return sum_chargefrac/((sumcharge/ndoms)*d_o)

def LikelihoodFunctor(data,domsUsed,vertexrad):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed

    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 18.949132224466762                                        # time parameter that has to be fit using simulations or data      
    vertexRad = vertexrad

    def likelihoodFunction(vtheta, vphi, theta, phi, t0):
	vertex = dataclasses.I3Position(vertexRad*np.sin(vtheta)*np.cos(vphi),vertexRad*np.sin(vtheta)*np.sin(vphi),vertexRad*np.cos(vtheta))
	vertex2 = dataclasses.I3Position(vertexRad*np.sin(theta)*np.cos(phi),vertexRad*np.sin(theta)*np.sin(phi),vertexRad*np.cos(theta))
	#direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))
        direction = dataclasses.I3Direction(vertex2.x-vertex.x,vertex2.y-vertex.y,vertex2.z-vertex.z)
        dark = 1.e-8

        N = 0.0;
        nsamples = 0.0
        sum_nloglike = 0.0
        for dom in pulse_series.keys() :
           domkey =  OMKey(dom.string, dom.om, 0) 
           d,dc,t,pdir = GetGeoTime(geo_doms[domkey].position,vertex,direction)
           pmtaccept = GetPMTAcceptance(pdir,dom.pmt);
           p_charge = np.exp(-d/tau)/max(dc,0.25)
           p_charge = p_charge*pmtaccept
           totalcharge = 0.0
           for pulse in pulse_series[dom] :
               charge = 1.0
               cpandel_out = cpandel(pulse.time - t0 - t ,d)
               if type(pulse_series) == type(dataclasses.I3RecoPulseSeriesMap()) :
                   charge = pulse.charge
               totalcharge += charge
               #sum_nloglike -= charge*np.log(cpandel_out+dark)
               sum_nloglike -= charge*np.log(cpandel_out*p_charge+dark)
               sum_nloglike -= min(0.0,pulse.time - t0 - t)
           N += totalcharge/(p_charge+dark)
           nsamples += 1.0

        #N /= nsamples
        #for dom in geo_doms.keys() :
        #   domkey =  OMKey(dom.string, dom.om, 0)
        #   if domkey in pulse_series.keys() :
        #        continue
        #   d,dc,t = GetGeoTime(geo_doms[dom].position,vertex,direction)
        #   p_charge = np.exp(-d/tau)/max(dc,0.25)
        #   sum_nloglike += N*p_charge 

        return sum_nloglike

    return likelihoodFunction


parser = argparse.ArgumentParser()                                              
parser.add_argument("-o", "--outfile",type = str, default="./test_output.root",help="")
parser.add_argument("-i", "--infile",type=str, default="./test_input.i3",help="")
#parser.add_argument("-g", "--gcdfile",type=str, default=os.getenv('PONESRCDIR')+"/GCD/PONE_Phase1.i3.gz",help="")
parser.add_argument("-g", "--gcdfile",type=str,default="",help="")  

args = parser.parse_args()    

fout = ROOT.TFile("MuonRecoRes.root","RECREATE")  
DirectionRes_linefit = ROOT.TH1F("DirectionRes_linefit","",101,-1.0,1.02)
DirectionRes_llhfit = ROOT.TH1F("DirectionRes_llhfit","",101,-1.0,1.02)
VertexRes_llhfit = ROOT.TH1F("VertexRes_llhfit","",100,0.0,200.0)
VertexRes_linefit = ROOT.TH1F("VertexRes_linefit","",100,0.0,200.0)
DirectionRes_linefit_vs_minr = ROOT.TH2F("DirectionRes_linefit_vs_minr","",110,0.,550.,100,-1.0,1.0)
DirectionRes_llhfit_vs_minr = ROOT.TH2F("DirectionRes_llhfit_vs_minr","",110,0.,550.,100,-1.0,1.0)
ThetaResolution_linefit = ROOT.TH1F("ThetaResolution_linefit","",100,-1.5,1.5)
ThetaResolution_llhfit = ROOT.TH1F("ThetaResolution_llh","",100,-1.5,1.5)
ThetaResolution_prefit = ROOT.TH1F("ThetaResolution_prefit","",100,-1.5,1.5)
ThetaResolution_vs_TotalCharge_linefit = ROOT.TH2F("ThetaResolution_vs_TotalCharge_linefit","",100,0,1000,100,-1.5,1.5)
ThetaResolution_vs_TotalCharge_llhfit = ROOT.TH2F("ThetaResolution_vs_TotalCharge_llhfit","",100,0,1000,100,-1.5,1.5)
ThetaResolution_vs_NDoms_linefit = ROOT.TH2F("ThetaResolution_vs_NDoms_linefit","",100,0,100,100,-1.5,1.5)
ThetaResolution_vs_NDoms_llhfit = ROOT.TH2F("ThetaResolution_vs_NDoms_llhfit","",100,0,100,100,-1.5,1.5)
PhiResolution_linefit = ROOT.TH1F("PhiResolution_linefit","",100,-1.5,1.5)
PhiResolution_llhfit = ROOT.TH1F("PhiResolution_llffit","",100,-1.5,1.5)
PhiResolution_prefit = ROOT.TH1F("PhiResolution_prefit","",100,-1.5,1.5)
TimeResolution_llhfit = ROOT.TH1F("TImeResolution_llhfit","",2000,-1000,1000)
llhfit_time = ROOT.TH1F("llhfit_time","",1000,4000,8000)
Muon_time = ROOT.TH1F("muon_time","",1000,4000,8000)
zenith_truth = ROOT.TH1F("zenith_truth","",100,0.0,3.14)
zenith_linefit = ROOT.TH1F("zenith_linefit","",100,0.0,3.14)
zenith_llhfit = ROOT.TH1F("zenith_llhfit","",100,0.0,3.14)
TimeResolution_llhfit_fmaxpe = ROOT.TH2F("TImeResolution_llhfit_fmaxpe","",100,0.,1.0,2000,-1000,1000)

#LikelihoodSpace_TrueVertex = ROOT.TH2F()

_dir = "/data/p-one/tmcelroy/muons/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if ('.i3.gz' in x and 'Reco_Corsika__PhotonProp' in x)]      
graphcount = 0
domsUsed = None
n = 1.34
ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
c = 0.299
c_n = c/ngroup                                     # light in water
theta_c = np.arccos(1./n) 
gcd_file = dataio.I3File("/home/users/tmcelroy/pone_offline/GCD/PONE_Phase1.i3.gz")
for frame in gcd_file:
    domsUsed = frame['I3Geometry'].omgeo

for infile in file_list :
	print(infile)
	infilei3 = dataio.I3File(os.path.join(_dir,infile))  

	for frame in infilei3:                                                          
                                                                                
		if not frame.Has('I3EventHeader') :                                          
			continue

		MMCTrackList = frame['MMCTrackList']
		Muon = MMCTrackList[0].GetI3Particle()
		t_offset = float(frame["TimeShiftedMCPEMap_toffset"].value)
		pulse_series = frame["SignificanHits"]
		mcpulses = frame["TimeShiftedMCPEMap"]
		likelihood = LikelihoodFunctor(pulse_series,domsUsed,550.)

		ndoms = 0
		totalcharge = 0.0
		fmaxPE = []
		for dom in pulse_series.keys() :
			ndoms += 1
			domtotalcharge = 0.0
                	for pulse in pulse_series[dom] :
                        	domtotalcharge += pulse.charge
			totalcharge +=domtotalcharge
			fmaxPE.append(domtotalcharge)
		
		for i in range(len(fmaxPE)) :
			fmaxPE[i] /= totalcharge

		muon_direction = dataclasses.I3Direction(Muon.dir.x,Muon.dir.y,Muon.dir.z)
		p_2 = Muon.pos.x**2.0+Muon.pos.y**2.0+Muon.pos.z**2.0
		pd = (Muon.pos.x*muon_direction.x+Muon.pos.y*muon_direction.y+Muon.pos.z*muon_direction.z)
		r_2 = 550.0**2.0

		L1 = -pd - np.sqrt(pd**2.0-p_2+r_2)
		L2 = -pd + np.sqrt(pd**2.0-p_2+r_2)
		muon_vertex = dataclasses.I3Position(Muon.pos.x+L1*muon_direction.x,Muon.pos.y+L1*muon_direction.y,Muon.pos.z+L1*muon_direction.z)
		muon_vertex2 = dataclasses.I3Position(Muon.pos.x+L2*muon_direction.x,Muon.pos.y+L2*muon_direction.y,Muon.pos.z+L2*muon_direction.z) 
		muon_time = Muon.time -t_offset +L1/0.299792458

		#print("vertex = "+str(Muon.pos.x)+","+str(Muon.pos.y)+","+str(Muon.pos.z))
		#print("direction = "+str(muon_direction.x)+","+str(muon_direction.y)+","+str(muon_direction.z))
		#print("p_r = "+str(np.sqrt(p_2)))
		#print("L = "+str(L))
		#print("t_offset = "+str(t_offset))
		#print("muon.time = "+str(Muon.time))
		#print("muontime = "+str(muon_time))

		#Closest Approach
		minl = -pd
		minpos = dataclasses.I3Position(Muon.pos.x+minl*muon_direction.x,Muon.pos.y+minl*muon_direction.y,Muon.pos.z+minl*muon_direction.z)
		minr = np.sqrt(minpos.x**2+minpos.y**2+minpos.z**2)

		#within detector
		min_cyl_r = np.sqrt(minpos.x**2+minpos.y**2)
		min_z = minpos.z

                linefit = frame['linefit']
                llhfit = frame['llhfit']
		prefit = frame['llhfit_prefit']

                p_2 = linefit.pos.x**2+linefit.pos.y**2+linefit.pos.z**2
                pd = (linefit.pos.x*linefit.dir.x+linefit.pos.y*linefit.dir.y+linefit.pos.z*linefit.dir.z)
                r_2 = 550.0**2

                L1 = -pd - np.sqrt(pd**2.0-p_2+r_2)
                L2 = -pd + np.sqrt(pd**2.0-p_2+r_2)

                linefit_vertex = dataclasses.I3Position(linefit.pos.x+L1*linefit.dir.x,linefit.pos.y+L1*linefit.dir.y,linefit.pos.z+L1*linefit.dir.z)
		linefit_vertex2 = dataclasses.I3Position(linefit.pos.x+L2*linefit.dir.x,linefit.pos.y+L2*linefit.dir.y,linefit.pos.z+L2*linefit.dir.z)
                linefit_dir = dataclasses.I3Position(linefit.dir.x,linefit.dir.y,linefit.dir.z)

		p_2 = llhfit.pos.x**2+llhfit.pos.y**2+llhfit.pos.z**2
                pd = (llhfit.pos.x*llhfit.dir.x+llhfit.pos.y*llhfit.dir.y+llhfit.pos.z*llhfit.dir.z)
                r_2 = 550.0**2

                L2 = max(-pd + np.sqrt(pd**2.0-p_2+r_2),-pd - np.sqrt(pd**2.0-p_2+r_2))
                llhfit_vertex2 = dataclasses.I3Position(llhfit.pos.x+L2*llhfit.dir.x,llhfit.pos.y+L2*llhfit.dir.y,llhfit.pos.z+L2*llhfit.dir.z)
                #print(np.sqrt(llhfit_vertex2.x**2.0+llhfit_vertex2.y**2.0+llhfit_vertex2.z**2.0))

		if totalcharge > 600 and graphcount<1:
		#if(totalcharge > 100 and totalcharge < 200 and 300 < np.abs(llhfit.time-muon_time) and graphcount<1):
		#	print("new graphs")
			graphcount += 1
			fout.cd()
			LikelihoodSpace_TrueVertex = ROOT.TH2F("Event_"+str(muon_vertex2.theta)+" "+str(muon_vertex2.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        LikelihoodSpace_TrueDirection = ROOT.TH2F("Event_"+str(muon_vertex.theta)+" "+str(muon_vertex.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        LikelihoodSpace_TrueVertandDir = ROOT.TH1F("Event_"+str(muon_time),"",8000,-2000.,2000.)
                        LikelihoodSpace_LineVertex = ROOT.TH2F("Eventline_"+str(llhfit_vertex2.theta)+" "+str(llhfit_vertex2.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        LikelihoodSpace_LineDirection = ROOT.TH2F("Eventline_"+str(llhfit.pos.theta)+" "+str(llhfit.pos.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        LikelihoodSpace_LineVertandDir = ROOT.TH1F("Eventline_"+str(llhfit.time),"",8000,-2000.,2000.)
                        QualitySpace_TrueVertex = ROOT.TH2F("QEvent_"+str(muon_vertex2.theta)+" "+str(muon_vertex2.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        QualitySpace_TrueDirection = ROOT.TH2F("QEvent_"+str(muon_vertex.theta)+" "+str(muon_vertex.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        QualitySpace_TrueVertandDir = ROOT.TH1F("QEvent_"+str(muon_time),"",8000,-2000.,2000.)
                        QualitySpace_LineVertex = ROOT.TH2F("QEventline_"+str(llhfit_vertex2.theta)+" "+str(llhfit_vertex2.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        QualitySpace_LineDirection = ROOT.TH2F("QEventline_"+str(llhfit.pos.theta)+" "+str(llhfit.pos.phi),"",120,0.,2.0*np.pi,60,0.,np.pi)
                        QualitySpace_LineVertandDir = ROOT.TH1F("QEventline_"+str(llhfit.time),"",8000,-2000.,2000.)
			time_10_20 = ROOT.TH1F("time_10_20_"+str(muon_time),"",3000,-3000.0,3000.0)
			time_50_60 = ROOT.TH1F("time_50_60_"+str(muon_time),"",3000,-3000.0,3000.0)
			time_80_100 = ROOT.TH1F("time_80_100_"+str(muon_time),"",3000,-3000.0,3000.0)
			pulsetime = ROOT.TH1F("pulsetime_"+str(muon_time),"",3000,-3000.0,3000.0)
			mcpulsetime  = ROOT.TH1F("mcpulsetime_"+str(muon_time),"",3000,-3000.0,3000.0)
			ChargeHits = ROOT.TGraph2D()
			mint_loglike = 999999999
			mint = 0
			for i in range(8000):
                                t0 = muon_time-2000.+0.5*float(i)
				likelihood_value = likelihood(muon_vertex.theta,muon_vertex.phi,muon_vertex2.theta,muon_vertex2.phi,t0)
				if likelihood_value < mint_loglike :
					mint_loglike = likelihood_value
					mint = t0
                                LikelihoodSpace_TrueVertandDir.SetBinContent(i+1,likelihood_value)

			for i in range(60) :
				for j in range(120) :
					#print(str(j)+","+str(i))
					phi = LikelihoodSpace_TrueVertex.GetXaxis().GetBinCenter(j+1)
					theta = LikelihoodSpace_TrueVertex.GetYaxis().GetBinCenter(i+1)
					_vertex = dataclasses.I3Position(550.*np.sin(theta)*np.cos(phi),550.*np.sin(theta)*np.sin(phi),550.*np.cos(theta))
					_direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))
					LikelihoodSpace_TrueVertex.SetBinContent(j+1,i+1,likelihood(muon_vertex.theta,muon_vertex.phi,theta, phi,muon_time))
					LikelihoodSpace_TrueDirection.SetBinContent(j+1,i+1,likelihood(theta,phi,muon_vertex2.theta,muon_vertex2.phi,muon_time))
					QualitySpace_TrueVertex.SetBinContent(j+1,i+1,timeresidualvariance(muon_vertex,_direction,pulse_series,domsUsed))
                                        QualitySpace_TrueDirection.SetBinContent(j+1,i+1,timeresidualvariance(_vertex,muon_direction,pulse_series,domsUsed))


                        mint_loglike = 999999999
                        for i in range(8000):
                                t0 = llhfit.time-2000.+0.5*float(i)
                                likelihood_value = likelihood(llhfit.pos.theta,llhfit.pos.phi,llhfit_vertex2.theta,llhfit_vertex2.phi,t0)
                                if likelihood_value < mint_loglike :
                                        mint_loglike = likelihood_value
                                        mint = t0
                                LikelihoodSpace_LineVertandDir.SetBinContent(i+1,likelihood_value)

                        for i in range(60) :
                                for j in range(120) :
                                        #print(str(j)+","+str(i))
                                        phi = LikelihoodSpace_TrueVertex.GetXaxis().GetBinCenter(j+1)
                                        theta = LikelihoodSpace_TrueVertex.GetYaxis().GetBinCenter(i+1)
					_vertex = dataclasses.I3Position(550.*np.sin(theta)*np.cos(phi),550.*np.sin(theta)*np.sin(phi),550.*np.cos(theta))
                                        _direction = dataclasses.I3Direction(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))
                                        LikelihoodSpace_LineVertex.SetBinContent(j+1,i+1,likelihood(llhfit.pos.theta,llhfit.pos.phi,theta, phi,llhfit.time))
                                        LikelihoodSpace_LineDirection.SetBinContent(j+1,i+1,likelihood(theta,phi,llhfit_vertex2.theta,llhfit_vertex2.phi,llhfit.time))
					QualitySpace_LineVertex.SetBinContent(j+1,i+1,timeresidualvariance(llhfit.pos,_direction,pulse_series,domsUsed))
                                        QualitySpace_LineDirection.SetBinContent(j+1,i+1,timeresidualvariance(_vertex,llhfit.dir,pulse_series,domsUsed))

			for i in range(1,3001) :

				time_10_20.SetBinContent(i,cpandel(time_10_20.GetXaxis().GetBinCenter(i),25))
				time_50_60.SetBinContent(i,cpandel(time_50_60.GetXaxis().GetBinCenter(i),55))
				time_80_100.SetBinContent(i,cpandel(time_80_100.GetXaxis().GetBinCenter(i),90))		
	
			count = 0
			for dom in pulse_series.keys() :
                        	for pulse in pulse_series[dom] :
                                	totalcharge += pulse.charge
                                        #d,dc,t = GetGeoTime(domsUsed[dom].position,muon_vertex,llhfit.dir)
                                        d,dc,t,pdir = GetGeoTime(domsUsed[dom].position,llhfit.pos,llhfit.dir)
					pulsetime.Fill(pulse.time-llhfit.time-t,pulse.charge)
				ChargeHits.SetPoint(count,domsUsed[dom].position.x,domsUsed[dom].position.y,domsUsed[dom].position.z);
				count +=1
			for dom in mcpulses.keys() :
				for pulse in mcpulses[dom] :
					domkey =  OMKey(dom.string, dom.om, 0)
                                        d,dc,t,pdir = GetGeoTime(domsUsed[domkey].position,muon_vertex,muon_direction)
					mcpulsetime.Fill(pulse.time-muon_time-t)
			LikelihoodSpace_TrueVertex.Write()
			LikelihoodSpace_TrueDirection.Write()
			LikelihoodSpace_TrueVertandDir.Write()
                        LikelihoodSpace_LineVertex.Write()
                        LikelihoodSpace_LineDirection.Write()
                        LikelihoodSpace_LineVertandDir.Write()
			QualitySpace_TrueVertex.Write()
                        QualitySpace_TrueDirection.Write()
                        QualitySpace_TrueVertandDir.Write()
                        QualitySpace_LineVertex.Write()
                        QualitySpace_LineDirection.Write()
                        QualitySpace_LineVertandDir.Write()
			time_10_20.Scale(totalcharge/time_10_20.Integral())
			time_10_20.Write()
			time_50_60.Scale(totalcharge/time_50_60.Integral())
			time_50_60.Write()
			time_80_100.Scale(totalcharge/time_80_100.Integral())
			time_80_100.Write()
			pulsetime.Write()
			mcpulsetime.Write()
			ChargeHits.Write("Event_"+str(totalcharge))
		
			#print("end graphs")

		llhfit_dir = dataclasses.I3Position(llhfit.dir.x,llhfit.dir.y,llhfit.dir.z)

		solidangle_linefit = muon_direction.x*linefit_dir.x+muon_direction.y*linefit_dir.y+muon_direction.z*linefit_dir.z
		solidangle_llhfit = muon_direction.x*llhfit_dir.x+muon_direction.y*llhfit_dir.y+muon_direction.z*llhfit_dir.z
		vertexdiff_linefit = np.sqrt((muon_vertex.x-linefit_vertex.x)**2.0+(muon_vertex.y-linefit_vertex.y)**2.0+(muon_vertex.z-linefit_vertex.z)**2.0)
		vertexdiff_llhfit = np.sqrt((muon_vertex.x-llhfit.pos.x)**2.0+(muon_vertex.y-llhfit.pos.y)**2.0+(muon_vertex.z-llhfit.pos.z)**2.0)
	
                if totalcharge > 100 :
		    DirectionRes_linefit.Fill(solidangle_linefit)
		    DirectionRes_llhfit.Fill(solidangle_llhfit)
		    VertexRes_llhfit.Fill(vertexdiff_llhfit)
		    VertexRes_linefit.Fill(vertexdiff_linefit)
		    DirectionRes_linefit_vs_minr.Fill(minr,solidangle_linefit)
		    DirectionRes_llhfit_vs_minr.Fill(minr,solidangle_llhfit)
		    ThetaResolution_linefit.Fill(linefit_dir.theta-muon_direction.theta)
		    ThetaResolution_llhfit.Fill(llhfit_dir.theta-muon_direction.theta)
		    ThetaResolution_prefit.Fill(prefit.dir.theta-muon_direction.theta)
		    PhiResolution_linefit.Fill((linefit_dir.phi-muon_direction.phi)*np.sin(muon_direction.theta))
		    PhiResolution_llhfit.Fill((llhfit_dir.phi-muon_direction.phi)*np.sin(muon_direction.theta))
                    PhiResolution_prefit.Fill((prefit.dir.phi-muon_direction.phi)*np.sin(muon_direction.theta))
		    TimeResolution_llhfit.Fill(llhfit.time-muon_time)
                    TimeResolution_llhfit_fmaxpe.Fill(max(fmaxPE),llhfit.time-muon_time)
                    llhfit_time.Fill(llhfit.time)
                    Muon_time.Fill(muon_time)
		    ThetaResolution_vs_TotalCharge_linefit.Fill(totalcharge,linefit_dir.theta-muon_direction.theta)
		    ThetaResolution_vs_TotalCharge_llhfit.Fill(totalcharge,llhfit_dir.theta-muon_direction.theta)
		    ThetaResolution_vs_NDoms_linefit.Fill(ndoms,linefit_dir.theta-muon_direction.theta)
		    ThetaResolution_vs_NDoms_llhfit.Fill(ndoms,llhfit_dir.theta-muon_direction.theta)
		    zenith_truth.Fill(muon_direction.theta)
		    zenith_linefit.Fill(linefit_dir.theta)
		    zenith_llhfit.Fill(llhfit_dir.theta)


fout.cd()
DirectionRes_linefit.Write()
DirectionRes_llhfit.Write()
VertexRes_llhfit.Write()
VertexRes_linefit.Write()
DirectionRes_linefit_vs_minr.Write()
DirectionRes_llhfit_vs_minr.Write()
ThetaResolution_linefit.Write()
ThetaResolution_llhfit.Write()
ThetaResolution_prefit.Write()
ThetaResolution_vs_TotalCharge_linefit.Write()
ThetaResolution_vs_TotalCharge_llhfit.Write()
ThetaResolution_vs_NDoms_linefit.Write()
ThetaResolution_vs_NDoms_llhfit.Write()
PhiResolution_linefit.Write()
PhiResolution_llhfit.Write()
PhiResolution_prefit.Write()
TimeResolution_llhfit.Write()
llhfit_time.Write() 
Muon_time.Write()
TimeResolution_llhfit_fmaxpe.Write()
zenith_truth.Write()
zenith_linefit.Write()
zenith_llhfit.Write()
fout.Close()










