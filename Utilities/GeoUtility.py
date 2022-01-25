import numpy as np
from icecube import dataio

"""!
Inputs:
    path to gcd file
Operation:
    Return cylinder radius and height based on input gcd file

"""

def get_geo_from_gcd(gcd):
    gcd_file = dataio.I3File(gcd)
    geometry = gcd_file.pop_frame()["I3Geometry"]
    doms = geometry.omgeo
    xPos = []
    yPos = []
    zPos = []
    for dom in doms.keys():
      x,y,z = doms[dom].position
      xPos.append(x)
      yPos.append(y)
      zPos.append(z)

    xPos = np.asarray(xPos)
    yPos = np.asarray(yPos)

    height = np.max(zPos)-np.min(zPos)
    rPos = np.sqrt(xPos*xPos+yPos*yPos)
    radius = np.max(rPos)
    gcd_file.close()
    return radius, height

