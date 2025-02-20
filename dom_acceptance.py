"""
Test script for comparing the old and new PMT acceptance code.
"""

from Utilities import DOMUtility
import numpy as np
from scipy.stats import rayleigh
import matplotlib.pyplot as plt
import os
import scipy.constants as const
from scipy.interpolate import CubicSpline, interp1d
from dataclasses import dataclass
from icecube import dataclasses, dataio, simclasses, icetray

dom_properties = DOMUtility.DOMProperties()


@dataclass
class Position:
    x: float
    y: float
    z: float

    @property
    def magnitude(self) -> float:
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)


@dataclass
class Photon:
    pos: Position
    wavelength: float


class POM:
    """
    Class characterizing the P-OM. COPY OF JAKUB's K40 code
    """

    def __init__(self, data_path):
        """
        Define P-OM properties and geometry
        """
        # ----------------------------------------------------------------------------
        # set pmt geometry
        # ----------------------------------------------------------------------------
        self.MODULE_RADIUS_M = 0.2159
        self.PMT_RADIUS_M = 0.055  # was 0.0635 before?
        # self.PMT_RADIUS_M = 0.0381

        PMT_ZENITHS = [32.5, 65.0, 115.0, 147.5]
        PMT_AZIMUTHS_TOP = [0.0, 90.0, 180.0, 270.0]
        PMT_AZIMUTHS_BOTTOM = [45.0, 135.0, 225.0, 315.0]

        zenith_list = np.array(sorted(4 * PMT_ZENITHS))
        azimuth_list = (
            PMT_AZIMUTHS_TOP
            + PMT_AZIMUTHS_BOTTOM
            + PMT_AZIMUTHS_BOTTOM
            + PMT_AZIMUTHS_TOP
        )

        x_coordinates = np.multiply(
            np.sin(np.deg2rad(zenith_list)), np.cos(np.deg2rad(azimuth_list))
        )
        y_coordinates = np.multiply(
            np.sin(np.deg2rad(zenith_list)), np.sin(np.deg2rad(azimuth_list))
        )
        z_coordinates = np.cos(np.deg2rad(zenith_list))

        self.PMT_MATRIX = np.array([x_coordinates, y_coordinates, z_coordinates]).T

        # ----------------------------------------------------------------------------
        # set efficiencies
        # ----------------------------------------------------------------------------

        self.COLLECTION_EFFICIENCY = 0.9
        self.TRANSIT_TIME_SPREAD = np.loadtxt(data_path + "tts.csv", delimiter=",")
        self.QUANTUM_EFFICIENCY = np.loadtxt(data_path + "qe.csv", delimiter=",")
        self.TRANSMITTANCE = np.loadtxt(data_path + "glass.csv", delimiter=",")
        self.ANGULAR_ACCEPTANCE_0 = np.loadtxt(data_path + "aa-0.csv", delimiter=",")

        energy_conversion = const.h * const.c / const.elementary_charge
        self.QE_FUNCTION = CubicSpline(
            np.flip(energy_conversion / self.QUANTUM_EFFICIENCY.T[0], axis=0),
            np.flip(self.QUANTUM_EFFICIENCY.T[1], axis=0),
            extrapolate=False,
        )
        self.T_FUNCTION = CubicSpline(
            self.TRANSMITTANCE.T[0] * 1e-9, self.TRANSMITTANCE.T[1], extrapolate=False
        )
        self.AA_FUNCTION = interp1d(
            self.ANGULAR_ACCEPTANCE_0.T[0],
            self.ANGULAR_ACCEPTANCE_0.T[1],
            bounds_error=False,
            fill_value=0,
        )

    def PMTHit(self, photon):
        """
        Function that return a list of PMTs on the
        P-OM that were hit by a given photon
        """

        photon_position = photon.pos

        # deterine the angle between all the pmt vectors and the photon vector
        pmt_photon_angles = self.PMT_MATRIX.dot(
            np.array([photon_position.x, photon_position.y, photon_position.z])
            / photon_position.magnitude
        )

        # check if the ditance at the module radius between the pmt vector
        # and the photon vector falls within the PMT size, if so, mark as a hit
        # return a list of indices that have been 'hit'
        pmt_vector_distances = self.MODULE_RADIUS_M * np.sin(
            np.arccos(pmt_photon_angles)
        )
        pmts_hit = np.where(
            np.logical_and(
                pmt_vector_distances < self.PMT_RADIUS_M, pmt_photon_angles >= 0
            )
        )[0]

        if len(pmts_hit) == 1:
            hit_angle = pmt_photon_angles[pmts_hit]
            hit_distance = pmt_vector_distances[pmts_hit]
            return [pmts_hit[0], hit_angle, hit_distance]
        elif len(pmts_hit) > 1:
            raise IndexError("two pmts are being hit at the same time!", pmts_hit)
        else:
            return [None, None, None]

    def CollectionEfficiency(self, photon_list):
        """
        Returns an array of collection efficiencies for
        a passed array of photons
        """
        return np.ones(len(photon_list)) * self.COLLECTION_EFFICIENCY

    def TransitTimeSpread(self, photon_list):
        """
        Returns an array of transit time spreads
        for a passed array of photons
        """
        return np.random.choice(
            self.TRANSIT_TIME_SPREAD[:, 0],
            len(photon_list),
            p=self.TRANSIT_TIME_SPREAD[:, 1],
        )

    def QuantumEfficiency(self, photon_list):
        """
        Returns an array of quantum efficiencies for
        a passed array of photons
        """
        wavelengths = np.array([p.wavelength for p in photon_list])
        efficiencies = np.nan_to_num(self.QE_FUNCTION(wavelengths))
        return efficiencies

    def Transmittance(self, photon_list):
        """
        Returns an array of transmittance efficiencies
        for a passed array of photons
        """
        wavelengths = np.array([p.wavelength for p in photon_list])
        efficiencies = np.nan_to_num(self.T_FUNCTION(wavelengths))
        return efficiencies

    def AngularAcceptanceSimple(self, pmt_info_list):
        """
        Returns an array of angular acceptances
        for a passed array of pmt information
        """
        distances = np.hstack(pmt_info_list[:, 2])
        efficiencies = np.nan_to_num(self.AA_FUNCTION(distances))
        return efficiencies


"""
This function is a copy of the GetPMT function from PONEDOMLauncher.py. 
Just for convenience.
"""


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


def k40_acceptance(module, photon_list):
    """Copied from https://github.com/pone-software/k40/blob/cf7804c3d25d6039bcdc06ca89629e0d7cf9c558/simulation/utilities/scripts/EventBuilder.py#L338"""

    number_photons = len(photon_list)
    pmt_hit_mask = np.zeros(number_photons, dtype=bool)
    pmt_list = np.empty(
        (number_photons, 3), dtype=object
    )  #### Need to make n x 3!!!!!!!!!

    for i in np.arange(number_photons):
        hit_pmt_info = module.PMTHit(photon_list[i])

        if hit_pmt_info[0] is not None:
            pmt_hit_mask[i] = True
            pmt_list[i] = hit_pmt_info

    filtered_photons = photon_list[pmt_hit_mask]
    filtered_pmts = pmt_list[pmt_hit_mask]

    if len(filtered_photons) == 0:
        return None

    photons = filtered_photons
    pmt_infos = np.vstack(filtered_pmts)

    # calculate efficiencies and remove photons with 0 efficency
    collection_efficiencies = module.CollectionEfficiency(photons)
    quantum_efficiencies = module.QuantumEfficiency(photons)
    transmittance_efficiencies = module.Transmittance(photons)
    angular_efficiencies = module.AngularAcceptanceSimple(pmt_infos)
    photon_efficiencies = (
        collection_efficiencies
        * quantum_efficiencies
        * transmittance_efficiencies
        * angular_efficiencies
    )

    return np.sum(photon_efficiencies)


def compare_acceptance_codes(n=1000000):

    # Generate directions uniformly on a sphere
    cos_zens = np.random.uniform(-1, 1, n)
    azis = np.random.uniform(0, 2 * np.pi, n)
    # Evaluate at 450nm
    wl = 450

    pmt_ixs = np.empty_like(cos_zens)
    for i, (cos_zen, azi) in enumerate(zip(cos_zens, azis)):

        theta = np.arccos(cos_zen)
        # Convert to cartesian coordinates
        photon_dir = [
            np.cos(azi) * np.sin(theta),
            np.sin(azi) * np.sin(theta),
            np.cos(theta),
        ]
        # Call GetPMT to get the PMT number, -1 means no hit
        pmt_ixs[i] = GetPMT(photon_dir, wl, 1)

    # Calculate acceptance
    accepted = np.sum(pmt_ixs >= 0)
    print("Old code: ", accepted / n)

    pmt_acc = DOMUtility.Geant4PMTAcceptance(
        os.getenv("PONESRCDIR") + "/data/pmt_acc.npz"
    )
    rel_positions = np.array(
        [
            [np.sin(np.arccos(cz)) * np.cos(a), np.sin(np.arccos(cz)) * np.sin(a), cz]
            for cz, a in zip(cos_zens, azis)
        ]
    )
    wavelengths = np.full(n, wl)

    pmt_ixs_new = pmt_acc.check_pmt_hit(
        rel_positions, wavelengths, np.ones_like(wavelengths)
    )

    accepted = np.sum(pmt_ixs_new > 0)
    print("New code: ", accepted / n)

    pom = POM("/home/chaack/repos/k40/simulation/utilities/data/")

    photon_list = []

    for rel_pos in rel_positions:
        photon = Photon(Position(*rel_pos), wl * 1e-9)
        photon_list.append(photon)

    accepted = k40_acceptance(pom, np.asarray(photon_list))

    print("K40 code", accepted / n)


def plot_acceptance_prob_map():

    # dom_properties.MakePMTAcceptancePlots(".")
    pmt_acc = np.load(os.getenv("PONESRCDIR") + "/data/ang_acceptance_16_all.npy")

    plt.figure()
    plt.imshow(
        np.sum(pmt_acc, axis=0).T / dom_properties.maxAngularAcceptance,
        extent=[0, 2 * np.pi, 0, np.pi],
        cmap=plt.cm.viridis,
        origin="lower",
        aspect="auto",
    )
    plt.xlabel("Phi (rad)")
    plt.ylabel("Theta (rad)")
    plt.colorbar(label="Bin Entry")
    plt.savefig("PMTAcceptance_sum.png")


def analyze_k40(file):
    hdl = dataio.I3File(file)
    n_hits = 0
    frames = 0
    while hdl.more():
        fr = hdl.pop_frame()

        if fr.Stop == icetray.I3Frame.DAQ:
            photons = fr["I3Photon_pmtsplit"]
            frames += 1

            times = sorted([p.time for _, ps in photons for p in ps])
            deltas = np.diff(times)

            if len(times) == 1:
                n_hits += 1
            elif len(times) > 1:
                print(times, sum(deltas > 20) + 1)
                n_hits += sum(deltas > 20) + 1
    print(frames)
    print(n_hits / (frames * 1e-5))
    hdl.close()


if __name__ == "__main__":
    # analyze_k40("clsim-output-1899_daq.i3.gz")
    compare_acceptance_codes()
    # plot_acceptance_prob_map()
