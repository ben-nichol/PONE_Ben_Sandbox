#!/usr/bin/env python3

import argparse
import h5py
import LeptonWeighter as LW
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate event weights using nuSQuIDS atmospheric flux."
    )
    parser.add_argument("--lic", required=True, help="LeptonInjector config file (mcgen.lic)")
    parser.add_argument("--events", required=True, help="Event HDF5 file (mcgen.h5)")
    parser.add_argument("--nusquids_flux", required=True, help="nuSQuIDS flux file (AtmFlux_output.h5)")
    parser.add_argument("--output", default="nusquids_weights.txt", help="Output text file for weights")
    parser.add_argument('--cross-section-location', default='/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/', help='Directory containing cross section files')
    return parser.parse_args()

def check_nusquids_support():
    if not hasattr(LW, 'nuSQUIDSAtmFlux'):
        print("ERROR: LeptonWeighter was not built with nuSQuIDS support.")
        sys.exit(1)

def main():
    args = parse_args()
    check_nusquids_support()

    # Load generator(s)
    generators = LW.MakeGeneratorsFromLICFile(args.lic)
    # Load nuSQuIDS flux
    nusquids_flux = LW.nuSQUIDSAtmFlux(args.nusquids_flux)
    
    # Load cross sections
    cross_section_location = args.cross_section_location
    xs = LW.CrossSectionFromSpline(
        cross_section_location + "/dsdxdy_nu_CC_iso.fits",
        cross_section_location + "/dsdxdy_nubar_CC_iso.fits",
        cross_section_location + "/dsdxdy_nu_NC_iso.fits",
        cross_section_location + "/dsdxdy_nubar_NC_iso.fits"
    )
      

    weighter = LW.Weighter(nusquids_flux, xs, generators)


    # Read events using h5py
    with h5py.File(args.events, "r") as h5file:
        events = h5file["EventProperties"]
        n_events = events.shape[0]
        with open(args.output, "w") as fout:
            fout.write("EventIdx\tEnergy_GeV\tZenith\tAzimuth\tWeight\n")
            for i in range(n_events):
                event = events[i]
                LWevent = LW.Event()
                # Field indices may need adjustment depending on your HDF5 structure
                LWevent.energy = event[5]
                LWevent.zenith = event[6]
                LWevent.azimuth = event[7]
                LWevent.interaction_x = event[8]
                LWevent.interaction_y = event[9]
                LWevent.final_state_particle_0 = LW.ParticleType(int(event[10]))
                LWevent.final_state_particle_1 = LW.ParticleType(int(event[11]))
                LWevent.primary_type = LW.ParticleType(int(event[12]))
                LWevent.radius = event[15] if len(event) > 15 else 0
                LWevent.total_column_depth = event[14] if len(event) > 14 else 0
                LWevent.x = 0
                LWevent.y = 0
                LWevent.z = event[16] if len(event) > 16 else 0

                weight = weighter(LWevent)
                fout.write(f"{i}\t{LWevent.energy:.6g}\t{LWevent.zenith:.6g}\t{LWevent.azimuth:.6g}\t{weight:.6e}\n")

    print(f"Done. Wrote weights for {n_events} events to {args.output}")

if __name__ == "__main__":
    main()
