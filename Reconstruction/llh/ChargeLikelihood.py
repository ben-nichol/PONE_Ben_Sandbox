import numpy as np
#from scipy.stats import poisson


# 1% probability that a PMT sees a dark pulse within event.
darkprob = 0.01
theta_cher = np.arccos(1./1.35)
alpha = np.sin(theta_cher)
tan_cher = np.tan(theta_cher)

def ComputeDOMAccept() :
	dom_acceptance = []
	nsamples = 100
	normalize = 1.0;
	for i in range(nsamples):
		angle = (np.pi*i)/float(nsamples)
		value = 0.0
		min_theta = angle
		max_theta = np.pi
		dtheta = (max_theta-min_theta)/float(nsamples)
		for j in range(nsamples) :
			theta = min_theta + dtheta*j
			if theta == 0.0 : 
				continue
			min_phi = 0.0
			max_phi = np.arccos(np.tan(angle)/np.tan(theta))
			if theta > np.pi/2.0 :
				max_phi = np.pi/2.0
			#print("max_phi = %f" % (max_phi))
			dphi = (max_phi-min_phi)/float(nsamples)
			for k in range(nsamples) :
				phi = min_phi + dphi*k
				value += -2.0*dphi*dtheta*np.sin(theta)*np.cos(theta+np.pi/2.0-angle)*np.cos(phi)
		if i == 0 :
			normalize = value
		value /= normalize
		dom_acceptance.append(value)
		#print("%f %f" % (angle,value))

dom_acceptance = ComputeDOMAccept()

def GetDOMEff(angle) :
	global dom_acceptance

	bin = angle*100./np.pi
	if bin > 99 :
		return 0.0
	if bin < 0 :
		return 0.0
	return dom_acceptance[bin]

#Hit probability energy independant
def HitProb(vert_x,vert_y,vert_z,theta,phi,dom_x,dom_y,dom_z) :
	global darkprob
	global theta_cher
	global alpha
	global tan_cher

	#Track is defined by it's sky position, thus everything is flipped.
	track = np.array([-np.cos(phi)*np.sin(theta),-np.sin(phi)*np.sin(theta),-np.cos(theta)])
	#Vertex
	vert = np.array([vert_x,vert_y,vert_z])
	# Orientation of Face of DOM, this is trivial for now but will be needed for mDOMs
	dom_orientation = np.array([0.,0.,-1.])
	#position of DOM
	dom_p = np.array([dom_x,dom_y,dom_z])

	#Vector from vertex to DOM
	vert_to_dom = vert-dom_p
	#Dot product with track to give distance to closest approach from vertex. 
	dot = np.dot(vert_to_dom,track)
	#Vector from point of cloases approach to DOM
	perp = vert_to_dom - dot*track
	#distance between closest approach and DOM
	d = np.sqrt(np.dot(perp,perp))

	#emission point of Cherenkov Photon
	emission_point = vert+(dot-d_phot)*track
	#vector from emission point to DOM
	emission_vec = dom_p-emission_point
	#Distance Photon travels
	d_phot = np.dot(emission_vec,emission_vec)
	#normalize emission vector to length 1
	emission_vec = emission_vec[0]/d_phot
	#Impact angle of photon on DOM.
	impact_theta = np.arcsin(np.dot(emission_vec,dom_orientation))

	#Assume DOM has a radius of 0.25m for now, tracks through DOM are not considered. This avoids zeros. 
	d_phot = max(d_phot,0.25)	

	#return the probability of seeing a photon from this angle at this distance including constant for darknoise.
	return max(1.0,GetDOMEff(impact_theta)/d_phot + darkprob)

def nLogLikelihood(pmt,charge,vert_x,vert_y,vert_z,theta,phi) :
	dom_x = pmt[:,0] 
    dom_y = pmt[:,1] 
    dom_z = pmt[:,2] 

    nLogLSum = 0.0
    for i in range(len(x)) : 
    	if charge[i] > 0.0 :
    		nLogLSum += -np.log(charge[i]*HitProb(vert_x,vert_y,vert_z,theta,phi,dom_x[i],dom_y[i],dom_z[i]))
    return nLogLSum

#Charge probability
def ChargeProb(vert_x,vert_y,vert_z,theta,phi,dom_x,dom_y,dom_z,N,charge) :
	global darkprob
	global theta_cher
	global alpha
	global tan_cher

	#Track is defined by it's sky position, thus everything is flipped.
	track = np.array([-np.cos(phi)*np.sin(theta),-np.sin(phi)*np.sin(theta),-np.cos(theta)])
	#Vertex
	vert = np.array([vert_x,vert_y,vert_z])
	# Orientation of Face of DOM, this is trivial for now but will be needed for mDOMs
	dom_orientation = np.array([0.,0.,-1.])
	#position of DOM
	dom_p = np.array([dom_x,dom_y,dom_z])

	#Vector from vertex to DOM
	vert_to_dom = vert-dom_p
	#Dot product with track to give distance to closest approach from vertex.
	dot = np.dot(vert_to_dom,track)
	#Vector from point of cloases approach to DOM
	perp = vert_to_dom - dot*track
	#distance between closest approach and DOM
	d = np.sqrt(np.dot(perp,perp))

	#emission point of Cherenkov Photon
	emission_point = vert+(dot-d_phot)*track
	#vector from emission point to DOM
	emission_vec = dom_p-emission_point
	#Distance Photon travels
	d_phot = np.dot(emission_vec,emission_vec)
	#normalize emission vector to length 1
	emission_vec = emission_vec[0]/d_phot
	#Impact angle of photon on DOM.
	impact_theta = np.arcsin(np.dot(emission_vec,dom_orientation))

 	#Assume DOM has a radius of 0.25m for now, tracks through DOM are not considered. This avoids zeros.
	d_phot = max(d_phot,0.25)

	#Get the number of expected photons
	mean_phot = N*GetDOMEff(impact_theta)/d_phot + darkprob

	# return the poisson probability for seeing Charge based on mean_phot, this assumes that SPE charge is close enough to 1.0
	# first order, will need to be better. 
	#return poisson.ppf(charge,mean_phot)
	
