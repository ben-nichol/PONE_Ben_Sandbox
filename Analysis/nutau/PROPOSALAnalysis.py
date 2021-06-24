# coding: utf-8
from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np

file = dataio.I3File('TestPropUpdated_config_combo_trunk.i3.gz')
#file = dataio.I3File('TestPropUpdated_config_cvmfs_combo.i3.gz')
qframes = []

while file.more():
    qframes.append(file.pop_daq())

energyRatioAll = ([])
energyRatio = ([])
AzDiff = ([])
ZenDiff = ([])
TRUE = 0
FALSE = 0
numberOfDaughterTaus = ([])
NuTauTrue = 0
NuTauFalse = 0
j = 1
for frame in qframes:
    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    tau = dataclasses.I3MCTree.first_child(mctree, primary[0].id)
    tauDirection = tau.dir
    tauDecayProd = dataclasses.I3MCTree.get_daughters(mctree, tau.id)

    daughterTaus = ([])
    daughterTauNeutrinos = ([])
    daughtersAll = ([])

    for i in range(0, len(tauDecayProd)):
        #print('FRAME NUMBER',  j)
        decayProdDirections = tauDecayProd[i].dir
        if tauDecayProd[i].type == 15 or tauDecayProd[i].type == -15:
            #print('Daughter Taus',  tauDecayProd[i].type)
            daughterTaus = np.append(daughterTaus, tauDecayProd[i].type)
        if tauDecayProd[i].type == 16 or tauDecayProd[i].type == -16:
            daughterTauNeutrinos = np.append(daughterTauNeutrinos, tauDecayProd[i].type)
        daughtersAll = np.append(daughtersAll, tauDecayProd[i].energy)
        AzD = abs(tauDirection.azimuth - decayProdDirections.azimuth)
        AzDiff = np.append(AzDiff, AzD)
        ZenD = abs(tauDirection.zenith - decayProdDirections.zenith)
        ZenDiff = np.append(ZenDiff, AzD)

    sumPAll = sum(daughtersAll)
    energyRatioAll = np.append(energyRatioAll, sumPAll/tau.energy)

    numberOfDaughterTaus = np.append(numberOfDaughterTaus, len(daughterTaus))

    if len(daughterTaus) > 0:
        TRUE = TRUE+1        #True, has daughter Taus
    else:
        #print(j)
        FALSE = FALSE+1

    if len(daughterTaus) < 2:
        daughters = ([])
        for i in range(0, len(tauDecayProd)):
            daughters = np.append(daughters, tauDecayProd[i].energy)
        sumP = sum(daughters)
        energyRatio = np.append(energyRatio, sumP/tau.energy)

    LeptConsFDaughter = ([])
    if len(daughterTauNeutrinos) == 0:
        NuTauTrue += 1       #True, Lepton Number is not conserved
    else:
        NuTauFalse += 1
    j = j + 1

#plot Daughter Taus
labels = ['Tau decays into Tau', 'Tau doesnt decay into Tau']
t = [TRUE]
f = [FALSE]

x = np.arange(len(labels))
width = 0.35  # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(x[0], t, width, label='True')
rects2 = ax.bar(x[1], f, width, label='False')

ax.set_ylabel('Number of Events')
ax.set_title('Checking Tau Decay')
ax.set_xticks(x)
ax.set_xticklabels(labels)

plt.savefig('DaugterTaus_trunk.pdf', dpi = 200)
plt.clf()

plt.hist(numberOfDaughterTaus, bins = 100, log=True)
plt.xlabel('number of daughter taus produced')
plt.ylabel('Count')
plt.title('Number of Daughter Taus Produced in Tau Decay')

plt.savefig('NumDaughterTaus_trunk.pdf', dpi = 200)
plt.clf()

#Lepton Conservation
labels = ['Tau doesnt decau into NuTau', 'Tau decays into NuTau']
t = [NuTauTrue]
f = [NuTauFalse]

x = np.arange(len(labels))
width = 0.35  # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(x[0], t, width, label='True')
rects2 = ax.bar(x[1], f, width, label='False')

ax.set_ylabel('Number of Events')
ax.set_title('Lepton Number Conservation')
ax.set_xticks(x)
ax.set_xticklabels(labels)

plt.savefig('LeptonNumCons_trunk.pdf', dpi = 200)
plt.clf()

#Direction
plt.hist(AzDiff, bins = 500, log=True)
plt.xlabel('Difference')
plt.ylabel('Count')
plt.title('Azimuth')

plt.savefig('Azimuth_trunk.pdf', dpi = 200)
plt.clf()

plt.hist(ZenDiff, bins = 500, log=True)
plt.xlabel('Difference')
plt.ylabel('Count')
plt.title('Zenith')

plt.savefig('Zenith_trunk.pdf', dpi = 200)
plt.clf()

#Energy
plt.hist(energyRatioAll, bins = 100, log=True)
plt.xlabel('Ratio')
plt.ylabel('Count')
plt.title('Energy(including all events)')

plt.savefig('EnergyAll_trunk.pdf', dpi = 200)
plt.clf()

plt.hist(energyRatio, bins = 100, log=True)
plt.xlabel('Ratio')
plt.ylabel('Count')
plt.title('Energy(excluding events with more than one daugter tau)')

plt.savefig('EnergyEx_trunk.pdf', dpi = 200)
plt.clf()
