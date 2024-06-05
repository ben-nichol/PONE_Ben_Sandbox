from Utilities import DOMUtility
import numpy as np
import numpy as np
from scipy.stats import rayleigh
import matplotlib.pyplot as plt
import os

dom_properties = DOMUtility.DOMProperties()


def GetPMT(photonDir, wl, weight):
        random = np.random.uniform()

        theta = 0.0
        phi = 0.0
        if len(photonDir) < 3:
            theta = photonDir[0]
            phi = photonDir[1]
        else:
            theta = np.arccos(photonDir[2])
            phi = np.arctan2(photonDir[1], photonDir[0])

        while phi < 0.0:
            phi += 2.0 * np.pi

        thetaBin = max(0, min(178, int(180.0 * theta / np.pi)))
        phiBin = max(0, min(358, int(180.0 * phi / np.pi)))
        pmtprobs = []
        for i in range(len(dom_properties.PMTacceptance)):
            pmtprobs.append(
                dom_properties.PMTacceptance[i][thetaBin][phiBin]
                / dom_properties.maxAngularAcceptance
            )

        totalprob = sum(pmtprobs)
        if random > totalprob:
            # print("angular kill photon")
            return -1

        random = np.random.uniform()
        i = 0
        sumprob = pmtprobs[0]
        while random > sumprob / totalprob and i < len(pmtprobs) - 1:
            i += 1
            sumprob += pmtprobs[i]

        qe = dom_properties.GetPMTQEnm(wl)
        # note: the probaility should only be weight * qe. It also includes maxAngularAcceptance
        #       here, since we already took this factor into account in `GetPMT` (which also rejects
        #       hits.)
        probability = weight * qe * dom_properties.maxAngularAcceptance

        # print("probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format( probability, weight, qe, dom_properties.maxAngularAcceptance, wl ))

        if probability > 1.0:
            print(
                "ERROR: probability={}. Weight={}, QE={}, maxAngularAcceptance={}, wavelength={}".format(
                    probability,
                    weight,
                    qe,
                    dom_properties.maxAngularAcceptance,
                    wl,
                )
            )
            raise RuntimeError(
                "The combined detection probability should never be > 1. You need to re-generate your I3Photons with a higher overall bias weight (setting `WavelengthAcceptance`)"
            )

        ### reject photon according to `probability`
        random = np.random.uniform()
        if random > probability:
            # print("kill photon")
            return -1

        # print("pmt = "+str(i+1))

        return i + 1



n = 100000
# Generate directions uniformly on a sphere
cos_zens = np.random.uniform(-1, 1, n)
azis = np.random.uniform(0, 2*np.pi, n)
# Evaluate at 410nm
wl = 450

pmt_ixs = np.empty_like(cos_zens)
for i, (cos_zen, azi) in enumerate(zip(cos_zens, azis)):

    theta = np.arccos(cos_zen)
    # Convert to cartesian coordinates
    photon_dir = [np.cos(azi) * np.sin(theta), np.sin(azi) * np.sin(theta), np.cos(theta)]
    # Call GetPMT to get the PMT number, 0 means no hit
    pmt_ixs[i] = GetPMT(photon_dir, wl, 1)

# Calculate acceptance
accepted = np.sum(pmt_ixs >= 0)
print("Old code: ", accepted / n)

pmt_acc = DOMUtility.Geant4PMTAcceptance(os.getenv("PONESRCDIR") +"/data/pmt_acc.npz")
rel_positions = np.array([[np.sin(np.arccos(cz)) * np.cos(a), np.sin(np.arccos(cz)) * np.sin(a), cz] for cz, a in zip(cos_zens, azis)])
wavelengths = np.full(n, wl)

pmt_ixs_new = pmt_acc.check_pmt_hit(rel_positions, wavelengths, np.ones_like(wavelengths))

accepted = np.sum(pmt_ixs_new > 0)
print("New code: ", accepted / n)

#dom_properties.MakePMTAcceptancePlots(".")
pmt_acc = np.load(os.getenv("PONESRCDIR") + "/data/ang_acceptance_16_all.npy")

plt.figure()
plt.imshow(np.sum(pmt_acc, axis=0).T / dom_properties.maxAngularAcceptance, extent=[0, 2*np.pi, 0, np.pi], cmap=plt.cm.viridis, origin="lower", aspect="auto")
plt.xlabel("Phi (rad)")
plt.ylabel("Theta (rad)")
plt.colorbar(label="Bin Entry")
plt.savefig("PMTAcceptance_sum.png")

"""
from icecube import dataio, dataclasses, simclasses, icetray

hdl = dataio.I3File("GenerateSingleMuons_100_photonprop.i3.gz")
ix_photons_map = {}
while hdl.more():

    fr = hdl.pop_frame()
    if fr.Stop == icetray.I3Frame.DAQ:
        photons = fr["I3Photons"]
        ix_photons_map[fr["I3EventHeader"].event_id] = sum([len(ps) for key, ps in photons])
hdl.close()

ix_photons_split_map = {}
hdl = dataio.I3File("GenerateSingleMuons_100_photonprop_daqSim_noise_ON.i3.gz")
while hdl.more():
    
    fr = hdl.pop_frame()
    if fr.Stop == icetray.I3Frame.DAQ:
        pmt_split_map = fr["I3Photons_pmtsplit"]
        print(fr["I3Photons_pmtsplit"])
        ix_photons_split_map[fr["I3EventHeader"].event_id] = sum([len(ps) for key, ps in pmt_split_map])

hdl.close()

for key in ix_photons_map.keys():
    if key in ix_photons_split_map:
        print(key, ix_photons_map[key], ix_photons_split_map[key])
"""
