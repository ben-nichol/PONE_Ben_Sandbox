import numpy as np

"""!
rand_dir
Inputs:
    initdir = x,y,z or theta,phi for the initial direction
    alpha = angel to shift initial direction by
Operation:
    Takes an initial direction and rotates it in random direction by angle alpha. 
    This is used for testing reconstructions. The output mimics the corrdinates of initdir.

"""
def rand_dir(initdir, alpha):
    theta = 0.0
    phi = 0.0
    if len(initdir) < 3 :
        theta = initdir[0]
        phi = initdir[1]
    else :
        theta = np.arccos(initdir[2])
        phi = np.arctan2(initdir[1],initdir[0])

    theta = 2*np.pi*rand.uniform(0,1)
    ct_dir = np.cos(theta)
    st_dir = np.sin(theta)
    cp_dir = np.cos(phi)
    sp_dir = np.sin(phi)
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    ct = np.cos(theta)
    st = np.sin(theta)
    dir_x = cp_dir*ct_dir*sa*ct + st_dir*ca*cp_dir - sp_dir*sa*st
    dir_y = sp_dir*ct_dir*sa*ct + st_dir*ca*sp_dir + cp_dir*sa*st
    dir_z = ca*ct_dir - sa*ct*st_dir
    if len(initdir)<3 :
        return np.array([np.arccos(dir_z),np.arctan(dir_y,dir_x)])
    return np.array([dir_x, dir_y, dir_z])


"""!
GetGeoTime(position,vert,direction)
Inputs:
    position = position of the DOM (x,y,z)
    vert = vertex point of track
    direction = direction of track
Operation:
    Computes the distance the cherenkov photon traveled through the water, 
    the distance of closest approach of the track to the DOM, 
    and the time since track was at vertex to the photon hitting the DOM.
"""
def GetGeoTime(position,vert,direction) :
    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water
    theta_c = np.arccos(1./n)
    x = position[0] - vert[0]
    y = position[1] - vert[1]
    z = position[2] - vert[2]
    dotprod = 0.0
    if len(direction) < 3 :
        dotprod = x*np.sin(direction[0])*np.cos(direction[1]) + y*np.sin(direction[0])*np.sin(direction[1]) + z*np.cos(direction[0])
    else :
        dotprod = x*direction[0] + y*direction[1] + z*direction[2]

    emission_point = []
    dc = max(0.25,x*x + y*y + z*z-dotprod*dotprod)
    ed = (dc/np.tan(theta_c))
    emission_point = [ed*direction[0],ed*direction[1],ed*direction[2]]
    emission_dir = [(position[0]-emission_point[0])/dc,(position[1]-emission_point[1])/dc,(position[2]-emission_point[2])/dc]
    theta = np.arccos(-emission_dir[2])
    phi = np.arccos(-emission_dir[0]/np.sqrt(1.-emission_dir[2]**2.0))
    dc = np.sqrt(dc)
    d = dc/np.sin(theta_c)
    t = d/c_n + dotprod/c - ed*c
    return d,dc,t, theta, phi

"""!
GetPhotonTravelTime(position,vert)
Inputs:
    position = I3Position for the DOM position
    vert = I3Position for vertex position.
Operation:
    Computes the distance traveled and travel time for the 
    photon assuming an isotropic emmission from the vertex.

"""

def GetPhotonTravelTime(position,vert):
    c = 0.299792458                                 # speed of light 
    n = 1.34
    ngroup = 1.35557                                # 1.33 is the refractive index of water at 20 degrees C
    c_n = c/ngroup                                     # light in water

    x = position[0] - vert[0]
    y = position[1] - vert[1]
    z = position[2] - vert[2]
    dc = np.sqrt(x*x + y*y + z*z)
    t = dc/c_n
    return dc,t
