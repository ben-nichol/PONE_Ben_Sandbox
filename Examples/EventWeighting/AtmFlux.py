import numpy as np
import nuSQuIDS as nsq

def main():
    units = nsq.Const()
    NUM_NEUTRINOS = 3
    E_MIN_GEV = 1.0
    E_MAX_GEV = 1.0e4

    FLUX_SPECTRUM_INDEX = -2.0  # E^-2 spectrum

    NEUTRINO_TYPE = nsq.NeutrinoType.neutrino
    INCLUDE_INTERACTIONS = False
  
    #neutrino direction?
    COS_THETA = -1.0  # upgoing
    THETA_RAD = np.arccos(COS_THETA)

    # Neutrino mixing angles (radians)
    THETA_12 = 0.563942
    THETA_13 = 0.154085
    THETA_23 = 0.785398
    # Mass squared differences (eV^2)
    DM2_21 = 7.65e-05
    DM2_31 = 0.00247
    
    H_MAX_KM = 500.0 # maximum propagation distance in km
    GSL_STEP = nsq.GSL_STEP_RK4 # GSL stepper type
    REL_ERROR = 1.0e-5
    ABS_ERROR = 1.0e-5
    N0 = 1.0e18  # normalization constant for initial flux
    OUTPUT_FILE = "AtmFlux_output.h5"
    
    # ==== Run NuSQuIDS ====
    E_min = E_MIN_GEV * units.GeV
    E_max = E_MAX_GEV * units.GeV
    E_nodes = nsq.logspace(E_min, E_max, 200)
    nus = nsq.nuSQUIDS(E_nodes, NUM_NEUTRINOS, NEUTRINO_TYPE, INCLUDE_INTERACTIONS)

    earth_atm = nsq.EarthAtm()
    track_atm = earth_atm.MakeTrack(THETA_RAD)
    nus.Set_Body(earth_atm)
    nus.Set_Track(track_atm)

    nus.Set_MixingAngle(0, 1, THETA_12)
    nus.Set_MixingAngle(0, 2, THETA_13)
    nus.Set_MixingAngle(1, 2, THETA_23)

    nus.Set_SquareMassDifference(1, DM2_21)
    nus.Set_SquareMassDifference(2, DM2_31)

    nus.Set_h_max(H_MAX_KM * units.km)
    nus.Set_GSL_step(GSL_STEP)
    nus.Set_rel_error(REL_ERROR)
    nus.Set_abs_error(ABS_ERROR)
    nus.Set_ProgressBar(True)

    E_range = nus.GetERange()
    inistate = np.zeros((len(E_range), NUM_NEUTRINOS))
    for i, E in enumerate(E_range):
        inistate[i, 1] = N0 * E**(FLUX_SPECTRUM_INDEX)  # muon neutrino flux

    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    print("\nEvolving neutrino state...")
    nus.EvolveState()

    print("\nWriting outputs...")
    nus.WriteStateHDF5(OUTPUT_FILE)


if __name__ == "__main__":
    main()