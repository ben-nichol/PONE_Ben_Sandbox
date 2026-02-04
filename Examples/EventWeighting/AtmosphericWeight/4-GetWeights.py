#!/usr/bin/env python3

import argparse
import h5py
import LeptonWeighter as LW
import sys
import os

# Parse command-line arguments for multiple LIC and HDF5 files
def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate event weights using nuSQuIDS atmospheric flux with multiple LIC and HDF5 files."
    )
    parser.add_argument('--lic', nargs='+', required=True, help='One or more LeptonInjector config files (mcgen.lic)')
    parser.add_argument('--events', nargs='+', required=True, help='One or more Event HDF5 files (mcgen.h5)')
    parser.add_argument('--nusquids_flux', required=True, help='nuSQuIDS flux file (AtmFlux_output.h5)')
    parser.add_argument('--output', default='nusquids_weights.txt', help='Output text file for weights')
    parser.add_argument('--cross-section-location', default='/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/', help='Directory containing cross section files')
    return parser.parse_args()

# Check for nuSQuIDS support in LeptonWeighter
def check_nusquids_support():
    if not hasattr(LW, 'nuSQUIDSAtmFlux'):
        print("ERROR: LeptonWeighter was not built with nuSQuIDS support.")
        sys.exit(1)

 # Utility function to convert event properties to a LeptonWeighter Event and get the weight
    # Adjust field indices as needed for your HDF5 structure
def get_weight(event, weighter):
    LWevent = LW.Event()
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
    return weighter(LWevent)

def main():
    # Parse arguments and check support
    args = parse_args()
    check_nusquids_support()

    # Load cross section splines
    cross_section_location = args.cross_section_location
    xs = LW.CrossSectionFromSpline(
        cross_section_location + "/dsdxdy_nu_CC_iso.fits",
        cross_section_location + "/dsdxdy_nubar_CC_iso.fits",
        cross_section_location + "/dsdxdy_nu_NC_iso.fits",
        cross_section_location + "/dsdxdy_nubar_NC_iso.fits"
    )

    # Load nuSQuIDS atmospheric flux
    nusquids_flux = LW.nuSQUIDSAtmFlux(args.nusquids_flux)

    # Load all generators from the provided LIC files
    # Each event will be weighted against all generators
    generators = []
    for lic_file in args.lic:
        generators.extend(LW.MakeGeneratorsFromLICFile(lic_file))

    # Build the weighter using all generators and the atmospheric flux
    weighter = LW.Weighter(nusquids_flux, xs, generators)

    # Loop over all HDF5 event files and calculate weights for each event
    with open(args.output, "w") as fout:
        fout.write("HDF5_File\tEventIdx\tEnergy_GeV\tZenith\tAzimuth\tWeight\n")
        for h5_file in args.events:
            with h5py.File(h5_file, "r") as h5file:
                events = h5file["EventProperties"]
                n_events = events.shape[0]
                for i in range(n_events):
                    event = events[i]
                    weight = get_weight(event, weighter)
                    fout.write(f"{h5_file}\t{i}\t{event[5]:.6g}\t{event[6]:.6g}\t{event[7]:.6g}\t{weight:.6e}\n")
    print(f"Done. Wrote weights for all events to {args.output}")

if __name__ == "__main__":
    main()
