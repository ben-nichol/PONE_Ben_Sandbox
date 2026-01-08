# EventWeighting Example Scripts

This directory contains example scripts for generating configuration files, writing event data tables, and calculating event weights for LeptonInjector simulations.

## 1-WriteGeneratorConfig.py

- Reads a file containing the S frame with LeptonInjector configuration (`muons.i3`).
- Serializes the LeptonInjector configuration from the S frame into a `.lic` file (`config.lic`).

## 2-WriteDataTables.py

- Reads the same input file (`muons.i3`).
- Writes the Q frame event data (`I3EventHeader`, `EventProperties`, `I3MCTree`) to an HDF5 file (`output.h5`) using the `I3SimHDFWriter` segment.

## 3-GetEventWeights.py

- Loads the generator configuration from `config.lic`.
- Loads cross-section data from the IceCube cross-section data repository.
- Defines a power-law flux and builds a weighting object.
- Loads event properties from `output.h5`.
- Calculates and prints the weight for each event using the loaded flux, cross-sections, and generator configuration.

---

**Usage:**

1. Run ../1-GenerateMuons.py to produce muons.i3. May need to copy here or repath
2. Run `2-WriteGeneratorConfig.py` to generate the LeptonInjector configuration file (`config.lic`).
3. Run `3-WriteDataTables.py` to write event data tables to `output.h5`.
4. Run `4-GetEventWeights.py` to calculate and print event weights using the generated files.

