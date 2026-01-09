# EventWeighting Example Scripts

This directory contains example scripts for generating configuration files, writing event data tables, and calculating event weights for LeptonInjector simulations.


## 2-WriteGeneratorConfig.py

- Reads a file containing the S frame with LeptonInjector configuration (`muons.i3`).
- Serializes the LeptonInjector configuration from the S frame into a `.lic` file (`config.lic`).
- Writes Q frame event data (`I3EventHeader`, `EventProperties`, `I3MCTree`) to an HDF5 file (`output.h5`) using the `I3SimHDFWriter` segment.

## 3-GetEventWeights.py

- Loads the generator configuration from `config.lic`.
- Loads cross-section data from the IceCube cross-section data repository.
- Defines a power-law flux and builds a weighting object.
- Loads event properties from `output.h5`.
- Calculates and prints the weight for each event using the loaded flux, cross-sections, and generator configuration.

**Usage:**

1. Run `../1-GenerateMuons.py` (outside this folder) to produce `muons.i3`. You may need to copy or repath this file.
2. Run `2-WriteGeneratorConfig.py` to generate the LeptonInjector configuration file (`config.lic`) and write event data tables to `output.h5`.
3. Run `3-GetEventWeights.py` to calculate and print event weights using the generated files.

## Open Questions
How are and should multiple S frames handled?
Is the weighting correct?

