import numpy as np
import nuSQuIDS as nsq

def main():

    #CONSTANT DEFINITIONS
    units = nsq.Const()
    NUM_NEUTRINOS = 3
    E_MIN_GEV = 1.0
    E_MAX_GEV = 1.0e4

    FLUX_SPECTRUM_INDEX = -2.0  # E^-2 spectrum

    NEUTRINO_TYPE = nsq.NeutrinoType.both
    INCLUDE_INTERACTIONS = False
  
    ## Neutrino mixing angles (radians)
    THETA_12 = 0.563942
    THETA_13 = 0.154085
    THETA_23 = 0.785398
    ## Mass squared differences (eV^2)
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
    E_nodes = nsq.logspace(E_min, E_max, 100)
    Z_nodes = nsq.linspace(-1.0, 1.0, 40)

    nus = nsq.nuSQUIDSAtm(Z_nodes, E_nodes, NUM_NEUTRINOS, NEUTRINO_TYPE, INCLUDE_INTERACTIONS)

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
    Z_range = nus.GetCosthRange()
    inistate = np.zeros((len(Z_range), len(E_range), 2, NUM_NEUTRINOS))
    for ci, cz in enumerate(Z_range):
        for ei, E in enumerate(E_range):
            for rho in range(2):
                for flv in range(NUM_NEUTRINOS):
                    inistate[ci, ei, rho, flv] = N0 * E**(FLUX_SPECTRUM_INDEX)

    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    print("\nEvolving neutrino state...")
    nus.EvolveState()

    print("\nWriting outputs...")
    nus.WriteStateHDF5(OUTPUT_FILE, True)  # True to overwrite if file exists


if __name__ == "__main__":
    main()
