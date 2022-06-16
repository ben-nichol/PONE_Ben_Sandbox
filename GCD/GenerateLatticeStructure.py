import numpy as np

def generateLaticeSpots(nstrings = 10):

    offset = np.pi*(1./6)
    anglediff = np.pi*(1./3)
    neighbourangles = [offset, anglediff+offset,2.0*anglediff+offset, 3.0*anglediff+offset, 4.0*anglediff+offset, 5.0*anglediff+offset]

    stringposx = [0.0]
    stringposy = [0.0]
    theta = [0.0]
    spacing = 50.0

    while len(stringposx) < nstrings :

        minradius = 1000000.0
        minradstring = 0
        minradstringneighbours = 10000
        for i in range(len(stringposx)) :
                nneighbours = 0
                rad = np.sqrt((stringposx[i])**2.0+(stringposy[i])**2.0)
                for j in range(len(stringposx)):
                        if i==j :
                                continue
                        dist = np.sqrt((stringposx[j]-stringposx[i])**2.0+(stringposy[j]-stringposy[i])**2.0)
                        if dist<spacing*1.2 :
                                nneighbours += 1
                if nneighbours < len(neighbourangles) and rad <= minradius :
                        if rad < minradius :
                                minradius = rad;
                                minradstring = i
                                minradstringneighbours = nneighbours
                        elif nneighbours < minradstringneighbours :
                                minradius = rad
                                minradstring = i
                                minradstringneighbours = nneighbours

        maxneighours = 0
        maxneighbourstring = 0
        for j in range(len(neighbourangles)) :
                newposx = stringposx[minradstring]+spacing*np.sin(neighbourangles[j])
                newposy = stringposy[minradstring]+spacing*np.cos(neighbourangles[j])

                nneighbours = 0
                overlap = False
                for k in range(len(stringposx)) :
                        dist = np.sqrt((newposx-stringposx[k])**2.0+(newposy-stringposy[k])**2.0)
                        if dist < spacing*0.8 :
                                nneighbours = 0
                                overlap = True
                        if dist<spacing*1.2 :
                                nneighbours += 1
                if nneighbours > maxneighours and not overlap:
                        maxneighours = nneighbours
                        maxneighbourstring = j
        stringposx.append(stringposx[minradstring]+spacing*np.sin(neighbourangles[maxneighbourstring]))
        stringposy.append(stringposy[minradstring]+spacing*np.cos(neighbourangles[maxneighbourstring]))
        theta.append(neighbourangles[maxneighbourstring])

    mean_x = sum(stringposx)/len(stringposx)
    mean_y = sum(stringposy)/len(stringposy)

    for i in range(len(stringposx)):
        stringposx[i] = (stringposx[i]-mean_x)/spacing
        stringposy[i] = (stringposy[i]-mean_y)/spacing

    return stringposx, stringposy, theta

