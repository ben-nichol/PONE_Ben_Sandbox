from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
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

def LikelihoodFunctor(data,domsUsed,vertexrad):
    # turn PMT locations and time hits into numpy arrays for easier numpy algebra
    pulse_series = data
    geo_doms = domsUsed
    time = []
    charge = []

    if(type(pulse_series) == 'icecube.dataclasses.I3RecoPulseSeriesMap') :
        for dom in pulse_series.keys() :
            for pulse in pulse_series[dom] :
                time.append(pulse.time)
                charge.append(pulse.charge)
    else :
        for dom in pulse_series.keys() :
            for pulse in pulse_series[dom] :
                time.append(pulse.time)
                charge.append(1.0)

    vx = 0.0
    vy = 0.0
    vz = 0.0
    v =np.array([0.0,0.0,0.0])

    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n)                       # Cherenkov angle in water in radians
    lambda_s = 120.                                 # scattering length of light for violet light
    lambda_a = 15.                                  # absorption length of light for violet light
    tau = 18.949132224466762                                        # time parameter that has to be fit using simulations or data      
    vertexRad = vertexrad
    # min time index for the first hit PMT
    min_index = np.argmin(time)

    # The computations from here on require we find the time and distance of closest approach, d_i,c and t_i,c
    def closestApproach():
        # Compute vec{r} - vec{x}
        dc = []
        tc = []
        lc = []

        for dom in pulse_series.keys() :
          for pulse in pulse_series[dom] :
            x = geo_doms[dom].position.x - vx
            y = geo_doms[dom].position.y - vy
            z = geo_doms[dom].position.z - vz
            # Compute (\vec{r} - vec{x}) dot \vec{v}
            dotprod = x*v[0] + y*v[1] + z*v[2]
            # Compute the final vector components
            #x = x - dotprod*v[0]
            #y = y - dotprod*v[1]
            #z = z - dotprod*v[2]
            # Compute t_i,c and d_i,c
            dc.append(np.sqrt(x*x + y*y + z*z-dotprod*dotprod))
            tc.append(dotprod/c)
            lc.append(dotprod)
        return dc, tc, lc

    # Given the linefit and other parameters 
    def computeResiduals(dc, tc, lc, t0):
        t = []
        d = []
        l = []
        for i in range(len(dc)) :
          # Now we find the time of the photon emission
          _tc = tc[i] - dc[i]/(np.tan(theta_c)*c)
          l.append(lc[i]-dc[i]/np.tan(theta_c))
          # The first component of the geometric time
          d.append(dc[i]/np.sin(theta_c))
          t_geo = d[-1]/c_n
          # Apply our offset time to find the "true" closest approach time array. Multiply by 1E9 to change to nanoseconds
          _tc += t0
          # The total geometric time
          t_geo = t_geo + _tc
        # Residual time is now the difference between the geometric time and the observed time. This won't work with just the Pandel Function
          t.append(time[i] - t_geo)
        return d, t, l

    # uses the prior defined functions to build a likelihood function that when given a track (linefit) will produce a negative loglikelihood value
    def likelihoodFunction(vtheta, vphi, theta, phi, t0):
        vx = vertexRad*np.sin(vtheta)*np.cos(vphi)
        vy = vertexRad*np.sin(vtheta)*np.sin(vphi)
        vz = vertexRad*np.cos(vtheta)
        v =np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        dc, tc, lc = closestApproach()
        d, t, l = computeResiduals(dc, tc, lc, t0)
        p_charge=[]
        dark = 1.e-8

        N = 0.0
        for dom in pulse_series.keys():
          totalcharge = 0.0
          for pulse in pulse_series[dom]:
             totalcharge += pulse.charge
          x = geo_doms[dom].position.x - vx
          y = geo_doms[dom].position.y - vy
          z = geo_doms[dom].position.z - vz
          # Compute (\vec{r} - vec{x}) dot \vec{v}
          dotprod = x*v[0] + y*v[1] + z*v[2]
          # Compute the final vector components
          #x = x - dotprod*v[0]
          #y = y - dotprod*v[1]
          #z = z - dotprod*v[2]
          # Compute t_i,c and d_i,c
          d_c = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
          l = dotprod - d_c/np.tan(theta_c)
          d_p = d_c/np.sin(theta_c)
          N += totalcharge/(np.exp(-d_p/tau)/max(d_c,0.25)+1000.*dark)


        for i in range(len(dc)) :
        #  if l[i]<start or l[i]>start+stop :
        #    p_charge.append(0.0)
        #  else :
        #    p_charge.append(np.exp(-d[i]/tau)/max(dc[i],0.5))
            p_charge.append(np.exp(-d[i]/tau)/max(dc[i],0.25))
        out = cpandel(t,d)
        sum_nloglike = 0.0
        for i in range(len(out)) :
            #sum_nloglike -= charge[i]*(np.log(N*out[i]+dark))
            sum_nloglike -= charge[i]*np.log(out[i]*p_charge[i]+dark/N)
            #sum_nloglike -= charge[i]*np.log(N*out[i]*p_charge[i])
        #likelihood of not seeing any light
        #assume chargetotal = sum(N*p_charge[i]) -> N = chargetotal/sum(p_charge[i])
        #Hit=[]
        #NotHit=[]
        #for dom in geo_doms.keys() :

          #x = geo_doms[dom].position.x - vx
          #y = geo_doms[dom].position.y - vy
          #z = geo_doms[dom].position.z - vz
          # Compute (\vec{r} - vec{x}) dot \vec{v}
          #dotprod = x*v[0] + y*v[1] + z*v[2]
          # Compute the final vector components
          #x = x - dotprod*v[0]
          #y = y - dotprod*v[1]
          #z = z - dotprod*v[2]
          # Compute t_i,c and d_i,c
          #d_c = np.sqrt(x*x + y*y + z*z)
          #l = dotprod - d_c/np.tan(theta_c) 
          #d_p = d_c/np.sin(theta_c)
          #if l < start or l > start+stop :
          #  if dom in pulse_series.keys() :
          #    Hit.append(1000.*dark)
          #  else :
          #    NotHit.append(1000.*dark)
          #else :
          #  if dom in pulse_series.keys() :
          #    Hit.append(np.exp(-d_p/tau)/max(d_c,0.5)+1000.*dark)
          #  else :
          #    NotHit.append(np.exp(-d_p/tau)/max(d_c,0.5)+1000.*dark)

        #N = sum(charge)/sum(Hit)
        #for p in NotHit :
        #  sum_nloglike += p*N 

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
ThetaResolution_vs_TotalCharge_linefit = ROOT.TH2F("ThetaResolution_vs_TotalCharge_linefit","",100,0,1000,100,-1.5,1.5)
ThetaResolution_vs_TotalCharge_llhfit = ROOT.TH2F("ThetaResolution_vs_TotalCharge_llhfit","",100,0,1000,100,-1.5,1.5)
ThetaResolution_vs_NDoms_linefit = ROOT.TH2F("ThetaResolution_vs_NDoms_linefit","",100,0,100,100,-1.5,1.5)
ThetaResolution_vs_NDoms_llhfit = ROOT.TH2F("ThetaResolution_vs_NDoms_llhfit","",100,0,100,100,-1.5,1.5)
#LikelihoodSpace_TrueVertex = ROOT.TH2F()

_dir = "/data/p-one/tmcelroy/muons/"#args.infile                                                              
file_list_aux = os.listdir(_dir)                                                
file_list = [x for x in file_list_aux if ('.i3.gz' in x and 'Reco_Corsika__PhotonProp' in x)]      
graphcount = 0
domsUsed = None
n = 1.34
ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
c_n = c/ngroup                                     # light in water
theta_c = np.arccos(1./n) 
gcd_file = dataio.I3File("/home/users/tmcelroy/pone_offline/GCD/PONE_Phase1.i3.gz")
for frame in gcd_file:
    domsUsed = frame['I3Geometry'].omgeo

for infile in file_list :
	infilei3 = dataio.I3File(os.path.join(_dir,infile))  

	for frame in infilei3:                                                          
                                                                                
		if not frame.Has('I3EventHeader') :                                          
			continue

		MMCTrackList = frame['MMCTrackList']
		Muon = MMCTrackList[0].GetI3Particle()
		t_offset = float(frame["TimeShiftedMCPEMap_toffset"].value)
		pulse_series = frame["SignificanHits"]
		likelihood = LikelihoodFunctor(pulse_series,domsUsed,550.)

		ndoms = 0
		totalcharge = 0.0

		for dom in pulse_series.keys() :
			ndoms += 1
			#print(dom)
			#print(str(domsUsed[dom].position.x)+","+str(domsUsed[dom].position.y)+","+str(domsUsed[dom].position.z))
                	for pulse in pulse_series[dom] :
                        	totalcharge += pulse.charge
				#print("		time = "+str(pulse.time)+" charge = "+str(pulse.charge))

		muon_direction = dataclasses.I3Direction(Muon.dir.x,Muon.dir.y,Muon.dir.z)
		p_2 = Muon.pos.x**2.0+Muon.pos.y**2.0+Muon.pos.z**2.0
		pd = (Muon.pos.x*muon_direction.x+Muon.pos.y*muon_direction.y+Muon.pos.z*muon_direction.z)
		r_2 = 550.0**2.0

		L = -pd - np.sqrt(pd**2.0-p_2+r_2)

		muon_vertex = dataclasses.I3Position(Muon.pos.x+L*muon_direction.x,Muon.pos.y+L*muon_direction.y,Muon.pos.z+L*muon_direction.z)
		muon_time = Muon.time -t_offset +L/0.299792458

		print("vertex = "+str(Muon.pos.x)+","+str(Muon.pos.y)+","+str(Muon.pos.z))
		print("direction = "+str(muon_direction.x)+","+str(muon_direction.y)+","+str(muon_direction.z))
		print("p_r = "+str(np.sqrt(p_2)))
		print("L = "+str(L))
		print("t_offset = "+str(t_offset))
		print("muon.time = "+str(Muon.time))
		print("muontime = "+str(muon_time))

		#Closest Approach
		minl = -pd
		minpos = dataclasses.I3Position(Muon.pos.x+minl*muon_direction.x,Muon.pos.y+minl*muon_direction.y,Muon.pos.z+minl*muon_direction.z)
		minr = np.sqrt(minpos.x**2+minpos.y**2+minpos.z**2)

		#within detector
		min_cyl_r = np.sqrt(minpos.x**2+minpos.y**2)
		min_z = minpos.z

		if(totalcharge > 600 and graphcount<5):
			print("new graphs")
			graphcount += 1
			fout.cd()
			LikelihoodSpace_TrueVertex = ROOT.TH2F("Event_"+str(muon_vertex.theta)+" "+str(muon_vertex.phi),"",60,0.,2.0*np.pi,30,0.,np.pi)
                        LikelihoodSpace_TrueDirection = ROOT.TH2F("Event_"+str(muon_direction.theta)+" "+str(muon_direction.phi),"",60,0.,2.0*np.pi,30,0.,np.pi)
                        LikelihoodSpace_TrueVertandDir = ROOT.TH1F("Event_"+str(muon_time),"",2000,-2000.,2000.)
			time_10_20 = ROOT.TH1F("time_10_20_"+str(muon_time),"",3000,6000.0,9000.0)
			time_50_60 = ROOT.TH1F("time_50_60_"+str(muon_time),"",3000,6000.0,9000.0)
			time_80_100 = ROOT.TH1F("time_80_100_"+str(muon_time),"",3000,6000.0,9000.0)
			ChargeHits = ROOT.TGraph2D()
			mint_loglike = 999999999
			mint = 0
			for i in range(2000):
				likelihood_value = likelihood(muon_vertex.theta,muon_vertex.phi,muon_direction.theta,muon_direction.phi,muon_time-2000.+2.*float(i))
				if likelihood_value < mint_loglike :
					mint_loglike = likelihood_value
					mint = muon_time-2000.+2.*float(i)
                                LikelihoodSpace_TrueVertandDir.SetBinContent(i+1,likelihood(muon_vertex.theta,muon_vertex.phi,muon_direction.theta,muon_direction.phi,muon_time-2000.+2.*float(i)))

			for i in range(30) :
				for j in range(60) :
					print(str(j)+","+str(i))
					phi = LikelihoodSpace_TrueVertex.GetXaxis().GetBinCenter(j+1)
					theta = LikelihoodSpace_TrueVertex.GetYaxis().GetBinCenter(i+1)
					LikelihoodSpace_TrueVertex.SetBinContent(j+1,i+1,likelihood(muon_vertex.theta,muon_vertex.phi,theta, phi,mint))
					LikelihoodSpace_TrueDirection.SetBinContent(j+1,i+1,likelihood(theta,phi,muon_direction.theta,muon_direction.phi,mint))
			
			count = 0
			for dom in pulse_series.keys() :
				x = geo_doms[dom].position.x - muon_vertex.x
            			y = geo_doms[dom].position.y - muon_vertex.y
            			z = geo_doms[dom].position.z - muon_vertex.z
            			# Compute (\vec{r} - vec{x}) dot \vec{v}
            			dotprod = x*muon_direction.x + y*muon_direction.y + z*muon_direction.z
            			# Compute the final vector components
            			d_c  = np.sqrt(x*x + y*y + z*z-dotprod*dotprod)
				d_p = d_c/np.sin(theta_c)		
                        	for pulse in pulse_series[dom] :
					if()
                                	totalcharge += pulse.charge
		
				ChargeHits.SetPoint(count,domsUsed[dom].position.x,domsUsed[dom].position.y,domsUsed[dom].position.z);
				count +=1
			LikelihoodSpace_TrueVertex.Write()
			LikelihoodSpace_TrueDirection.Write()
			LikelihoodSpace_TrueVertandDir.Write()
			ChargeHits.Write("Event_"+str(totalcharge))
		
			print("end graphs")
		#if np.abs(min_z) > 500. or min_cyl_r > 200. :
		#	continue

		linefit = frame['linefit']
		llhfit = frame['llhfit']

		p_2 = linefit.pos.x**2+linefit.pos.y**2+linefit.pos.z**2
                pd = (linefit.pos.x*linefit.dir.x+linefit.pos.y*linefit.dir.y+linefit.pos.z*linefit.dir.z)
                r_2 = 550.0**2

                L = -pd - np.sqrt(pd**2.0-p_2+r_2)

                linefit_vertex = dataclasses.I3Position(linefit.pos.x+L*linefit.dir.x,linefit.pos.y+L*linefit.dir.y,linefit.pos.z+L*linefit.dir.z)
		linefit_dir = dataclasses.I3Position(linefit.dir.x,linefit.dir.y,linefit.dir.z)
		llhfit_dir = dataclasses.I3Position(llhfit.dir.x,llhfit.dir.y,llhfit.dir.z)

		solidangle_linefit = muon_direction.x*linefit_dir.x+muon_direction.y*linefit_dir.y+muon_direction.z*linefit_dir.z
		solidangle_llhfit = muon_direction.x*llhfit_dir.x+muon_direction.y*llhfit_dir.y+muon_direction.z*llhfit_dir.z
		vertexdiff_linefit = np.sqrt((muon_vertex.x-linefit_vertex.x)**2.0+(muon_vertex.y-linefit_vertex.y)**2.0+(muon_vertex.z-linefit_vertex.z)**2.0)
		vertexdiff_llhfit = np.sqrt((muon_vertex.x-llhfit.pos.x)**2.0+(muon_vertex.y-llhfit.pos.y)**2.0+(muon_vertex.z-llhfit.pos.z)**2.0)

		DirectionRes_linefit.Fill(solidangle_linefit)
		DirectionRes_llhfit.Fill(solidangle_llhfit)
		VertexRes_llhfit.Fill(vertexdiff_llhfit)
		VertexRes_linefit.Fill(vertexdiff_linefit)
		DirectionRes_linefit_vs_minr.Fill(minr,solidangle_linefit)
		DirectionRes_llhfit_vs_minr.Fill(minr,solidangle_llhfit)
		ThetaResolution_linefit.Fill(linefit_dir.theta-muon_direction.theta)
		ThetaResolution_llhfit.Fill(llhfit_dir.theta-muon_direction.theta)
		ThetaResolution_vs_TotalCharge_linefit.Fill(totalcharge,linefit_dir.theta-muon_direction.theta)
		ThetaResolution_vs_TotalCharge_llhfit.Fill(totalcharge,llhfit_dir.theta-muon_direction.theta)
		ThetaResolution_vs_NDoms_linefit.Fill(ndoms,linefit_dir.theta-muon_direction.theta)
		ThetaResolution_vs_NDoms_llhfit.Fill(ndoms,llhfit_dir.theta-muon_direction.theta)


fout.cd()
DirectionRes_linefit.Write()
DirectionRes_llhfit.Write()
VertexRes_llhfit.Write()
VertexRes_linefit.Write()
DirectionRes_linefit_vs_minr.Write()
DirectionRes_llhfit_vs_minr.Write()
ThetaResolution_linefit.Write()
ThetaResolution_llhfit.Write()
ThetaResolution_vs_TotalCharge_linefit.Write()
ThetaResolution_vs_TotalCharge_llhfit.Write()
ThetaResolution_vs_NDoms_linefit.Write()
ThetaResolution_vs_NDoms_llhfit.Write()
fout.Close()










