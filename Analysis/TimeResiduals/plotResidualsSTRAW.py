#check Residuals

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib, csv
from statistics import mean
from icecube.phys_services import I3Calculator

def timeRes(file, wavelength):
    pfile = dataio.I3File(str(file))
    x, y, z, wlen, time = ([]), ([]), ([]), ([]), ([])
    weights = ([])
    frame = pfile.pop_frame()

    for pframe in pfile:
        photonMap = pframe['I3Photons']
        for modKey,photons in photonMap:
            for photon in photons:
                wlen = np.append(wlen, photon.wavelength)
                #if photon.wavelength / I3Units.nanometer >= 360 and photon.wavelength / I3Units.nanometer <= 370:
                if photon.wavelength / I3Units.nanometer >= wavelength-5 and photon.wavelength / I3Units.nanometer <= wavelength+5:
                    photonPos = photon.pos
                    x = np.append(x, photonPos.x)
                    y = np.append(y, photonPos.y)
                    z = np.append(z, photonPos.z)
                    time = np.append(time, photon.time)
                    weights = np.append(weights, photon.weight)

    print("length of data - ", len(x))
    distance = np.sqrt(x**2 + y**2 + (z - 107.66)**2)
    boolean1 = [(distance > 51) & (distance < 55)]
    boolean2 = [(distance > 68) & (distance < 72)]
    boolean3 = [(distance > 86) & (distance < 90)]
    #boolean = [(distance > 54.2) & (distance < 54.6)]
    #boolean = [(distance > 70.3) & (distance < 70.7)]
    #boolean = [(distance > 39.8) & (distance < 40.2)]
    #boolean = [(distance > 87.9) & (distance < 88.3)]
    selectDist1 = distance[boolean1]
    selectDist2 = distance[boolean2]
    selectDist3 = distance[boolean3]
    W1 = weights[boolean1]
    W2 = weights[boolean2]
    W3 = weights[boolean3]
    #print('Distance Length - ', len(selectDist1))
    selectTime1 = time[boolean1]
    selectTime2 = time[boolean2]
    selectTime3 = time[boolean3]
    velocity = (299792458/1.33)
    actualTime1 = (selectDist1/velocity)*1e9
    actualTime2 = (selectDist2/velocity)*1e9
    actualTime3 = (selectDist3/velocity)*1e9
    #actualTime = (88.1/velocity)*1e9
    tRes1 = selectTime1 - actualTime1
    tRes2 = selectTime2 - actualTime2
    tRes3 = selectTime3 - actualTime3
    print ('tRes Length', len(tRes1))
    #print("Tres- ", tRes)
    smear1 = np.random.normal(1, 3)
    smear2 = np.random.normal(1, 3)
    smear3 = np.random.normal(1, 3)
    #smear = 0

    return tRes1+smear1, tRes2+smear2, tRes3+smear3, W1, W2, W3, wlen

timeResidsUV1,timeResidsUV2,timeResidsUV3  = ([]),([]),([])
weights1, weights2, weights3 = ([]), ([]), ([])
wavelengthUV = ([])

for i in range(0, 20):
    #if i == 2 or i == 4 or i == 7 or i==8 or i == 15:
        #continue
    #infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_1_sDOM_90.i3.gz'
    #infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_1_right.i3.gz'
    infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_1_20000_2e5_spread_STRAWGeo.i3.gz'
    #infile = '/home/users/akatil/P-ONE/medium/SCPhotonsSTRAW.i3'
    tRes1, tRes2, tRes3, W1, W2, W3, wlen = timeRes(infile,365)
    #timeResidsUV1 = np.append(timeResids1, tRes2)
    #timeResidsUV2 = np.append(timeResids2, tRes2)
    #timeResidsUV3 = np.append(timeResids3, tRes3)
    timeResidsUV1 = np.append(timeResidsUV1, tRes2)
    timeResidsUV2 = np.append(timeResidsUV2, tRes2)
    timeResidsUV3 = np.append(timeResidsUV3, tRes3)
    weightsUV1 = np.append(weights1, W1)
    weightsUV2 = np.append(weights2, W2)
    weightsUV3 = np.append(weights3, W3)
    wavelengthUV = np.append(wavelengthUV, wlen)

timeResidsVio1,timeResidsVio2,timeResidsVio3  = ([]),([]),([])
weights1, weights2, weights3 = ([]), ([]), ([])
wavelengthV = ([])

for i in range(0, 20):
    #if i == 2 or i == 4 or i == 7 or i==8 or i == 15:
        #continue
    #infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_2_sDOM_90.i3.gz'
    #infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_2_right.i3.gz'
    infile = '/data/p-one/akatil/timeResiduals/generatePhotons_' + str(i) + '_2_20000_2e5_spread_STRAWGeo.i3.gz'
    print (i)
    #infile = '/home/users/akatil/P-ONE/medium/SCPhotonsSTRAW.i3'
    tRes1, tRes2, tRes3, W1, W2, W3, wlen = timeRes(infile,405)
    #timeResidsVio1 = np.append(timeResids1, tRes1)
    #timeResidsVio2 = np.append(timeResids2, tRes2)
    #timeResidsVio3 = np.append(timeResids3, tRes3)
    timeResidsVio1 = np.append(timeResidsVio1, tRes1)
    timeResidsVio2 = np.append(timeResidsVio2, tRes2)
    timeResidsVio3 = np.append(timeResidsVio3, tRes3)
    weightsVio1 = np.append(weights1, W1)
    weightsVio2 = np.append(weights2, W2)
    weightsVio3 = np.append(weights3, W3)
    wavelengthV = np.append(wavelengthV, wlen)

#bins = np.linspace(-10,200,111)
#plt.figure()
#plt.hist(timeResids, bins = bins,  histtype = 'step', log = True, density = True, label = 'Simulation')
#plt.step(time, hits, where = 'post', label = 'STRAW data')
#plt.xlabel("Time Residual (ns)")
#plt.ylabel("Normalized Count")
#plt.title("Time Residual Distribution")
#plt.legend()
#plt.show()

#plt.savefig('timeResidualsSTRAW_violet_Direct_200ns.pdf', dpi = 200)

#combinedUV1 = np.vstack((timeResidsUV1, weightsUV1)).T
#combinedUV2 = np.vstack((timeResidsUV2, weightsUV2)).T
#combinedUV3 = np.vstack((timeResidsUV3, weightsUV3)).T
#combinedV1 = np.vstack((timeResidsVio1, weightsVio1)).T
#combinedV2 = np.vstack((timeResidsVio2, weightsVio2)).T
#combinedV3 = np.vstack((timeResidsVio3, weightsVio3)).T

np.savetxt('/home/users/akatil/P-ONE/tResiduals/vioDOM_0_54_20000_2e5_spread_STRAWGeo.csv', timeResidsVio1, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/vioDOM_0_70_20000_2e5_spread_STRAWGeo.csv', timeResidsVio2, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/vioDOM_0_88_20000_2e5_spread_STRAWGeo.csv', timeResidsVio3, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/uvDOM_0_54_20000_2e5_spread_STRAWGeo.csv', timeResidsUV1, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/uvDOM_0_70_20000_2e5_spread_STRAWGeo.csv', timeResidsUV2, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/uvDOM_0_88_20000_2e5_spread_STRAWGeo.csv', timeResidsUV3, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/uvWlen_0_20000_2e5_spread_STRAWGeo.csv', wavelengthUV, delimiter=',')
np.savetxt('/home/users/akatil/P-ONE/tResiduals/vioWlen_0_20000_2e5_spread_STRAWGeo.csv', wavelengthV, delimiter=',')

#np.savetxt('new_viosDOM_90_54_5000.csv', timeResidsVio1, delimiter=',')
#np.savetxt('new_viosDOM_90_70_5000.csv', timeResidsVio2, delimiter=',')
#np.savetxt('new_viosDOM_90_88_5000.csv', timeResidsVio3, delimiter=',')
#np.savetxt('new_uvsDOM_90_54_5000.csv', timeResidsUV1, delimiter=',')
#np.savetxt('new_uvsDOM_90_70_5000.csv', timeResidsUV2, delimiter=',')
#np.savetxt('new_uvsDOM_90_88_5000.csv', timeResidsUV3, delimiter=',')
