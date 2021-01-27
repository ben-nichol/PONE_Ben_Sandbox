import numpy as np
#from scipy.stats import poisson


# 1% probability that a PMT sees a dark pulse within event.
darkprob = 0.01
theta_cher = np.arccos(1./1.35)
alpha = np.sin(theta_cher)
tan_cher = np.tan(theta_cher)

def dotprod(vec1, vec2) :
  return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]

def scalevec(amp,vec) :
  vec[0] *= amp
  vec[1] *= amp
  vec[2] *= amp

  return vec

def addvec(vec1,vec2,scale) :

  return [vec1[0]+scale*vec2[0],vec1[1]+scale*vec2[1],vec1[2]+scale*vec2[2]]

def ComputeDOMAccept() :
  dom_acceptance = list()
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
        value += -2.0*dphi*dtheta*np.sin(theta)*np.cos(theta+np.pi/2.0-angle) #*np.cos(phi)
    if i == 0 :
      normalize = value
    value /= normalize
    dom_acceptance.append(value)
    #print("%f %f" % (angle,value))
    sum_all = sum(dom_acceptance)
    for i in range(len(dom_acceptance)) :
      dom_acceptance[i] = dom_acceptance[i]/sum_all
  return dom_acceptance

dom_acceptance = ComputeDOMAccept()

def GetDOMEff(angle) :
  global dom_acceptance

  #mdom uniform coverage
  return 1.0

  if len(dom_acceptance) == 0 :
    dom_acceptance = ComputeDOMAccept()

  binnum = int(angle*100./np.pi)
  if binnum > 99 :
    return 0.0
  if binnum < 0 :
    return 0.0
  return dom_acceptance[binnum]

#Hit probability energy independant
def HitProb(vert_x,vert_y,vert_z,theta,phi,dom_x,dom_y,dom_z) :
  global darkprob
  global theta_cher
  global alpha
  global tan_cher

	#Track is defined by it's sky position, thus everything is flipped.
  track = [np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]
	#Vertex
  vert = [vert_x,vert_y,vert_z]
	# Orientation of Face of DOM, this is trivial for now but will be needed for mDOMs
  dom_orientation = np.array([0.,0.,1.])
	#position of DOM
  dom_p = [dom_x,dom_y,dom_z]

	#Vector from vertex to DOM
  vert_to_dom = addvec(vert,dom_p,-1.)
	#Dot product with track to give distance to closest approach from vertex. 
  dot = dotprod(vert_to_dom,track)
	#Vector from point of closest approach to DOM
  perp = addvec(vert_to_dom,track,-dot)
	#distance between closest approach and DOM
  d = np.sqrt(dotprod(perp,perp))
  
	#emission point of Cherenkov Photon
  emission_point = addvec(vert,track,(dot-d/np.tan(theta_cher)))
	#vector from emission point to DOM
  emission_vec = addvec(perp,track,d/np.tan(theta_cher))
	#Distance Photon travels
  d_phot = np.sqrt(dotprod(emission_vec,emission_vec))
	#normalize emission vector to length 1
  emission_vec = scalevec(1./d_phot,emission_vec)
	#Impact angle of photon on DOM.
 
  impact_theta = np.arcsin(dotprod(emission_vec,dom_orientation))

	#Assume DOM has a radius of 0.25m for now, tracks through DOM are not considered. This avoids zeros. 
  d_phot = max(d_phot,0.25)	

	#return the probability of seeing a photon from this angle at this distance including constant for darknoise.

  #print("impact = %f DOMAcc = %f d_phot=%d %f" % (impact_theta,GetDOMEff(impact_theta),d_phot,min(1.0,GetDOMEff(impact_theta)/d_phot+darkprob)))
  #return min(1.0,np.exp(-d_phot/30.)*GetDOMEff(impact_theta)/d_phot + darkprob)
  return GetDOMEff(impact_theta)/d_phot

def Likelihood(pmt,charge,vert_x,vert_y,vert_z,theta,phi) :
  dom_x = pmt[:,0] 
  dom_y = pmt[:,1] 
  dom_z = pmt[:,2] 

  prob = []
  for i in range(len(dom_x)) : 
       prob.append(HitProb(vert_x,vert_y,vert_z,theta,phi,dom_x[i],dom_y[i],dom_z[i]))
  return prob

#Charge probability
def ChargeProb(vert_x,vert_y,vert_z,theta,phi,dom_x,dom_y,dom_z,N,charge) :
	global darkprob
	global theta_cher
	global alpha
	global tan_cher

	#Track is defined by it's sky position, thus everything is flipped.
	track = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)])
	#Vertex
	vert = np.array([vert_x,vert_y,vert_z])
	# Orientation of Face of DOM, this is trivial for now but will be needed for mDOMs
	dom_orientation = np.array([0.,0.,1.])
	#position of DOM
	dom_p = np.array([dom_x,dom_y,dom_z])

	#Vector from vertex to DOM
	vert_to_dom = vert-dom_p
	#Dot product with track to give distance to closest approach from vertex.
	dot = dotprod(vert_to_dom,track)
	#Vector from point of cloases approach to DOM
	perp = vert_to_dom - scalevec(dot,track)
	#distance between closest approach and DOM
	d = np.sqrt(dotprod(perp,perp))

	#emission point of Cherenkov Photon
	emission_point = vert+scalevec((dot-d/np.tan(theta_cher)),track)
	#vector from emission point to DOM
	emission_vec = dom_p-emission_point
	#Distance Photon travels
	d_phot = np.sqrt(dotprod(emission_vec,emission_vec))
	#normalize emission vector to length 1
	emission_vec = scalevec(1./d_phot,emission_vec)
	#Impact angle of photon on DOM.
	impact_theta = np.arcsin(np.dot(emission_vec,dom_orientation))

 	#Assume DOM has a radius of 0.25m for now, tracks through DOM are not considered. This avoids zeros.
	d_phot = max(d_phot,0.25)

	#Get the number of expected photons
	mean_phot = N*GetDOMEff(impact_theta)/d_phot + darkprob

	# return the poisson probability for seeing Charge based on mean_phot, this assumes that SPE charge is close enough to 1.0
	# first order, will need to be better. 
	#return poisson.ppf(charge,mean_phot)
	
