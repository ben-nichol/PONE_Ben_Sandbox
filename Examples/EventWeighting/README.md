
# EventWeighting Example Scripts

This directory contains example scripts for generating configuration files, writing event data tables, and calculating event weights for LeptonInjector simulations. There are now two main subfolders:

## Structure
- `BasicWeight/`
	- `1-GenSimpleMC.py`: Generates a simple Monte Carlo sample.
	- `2-WriteGeneratorConfig.py`: Serializes LeptonInjector configuration and writes event tables.
	- `3-GetEventWeights.py`: Calculates event weights using a basic power-law flux.
	- `config.lic`, `muons.i3`, `output.h5`: Configuration, input, and output files.

- `AtmosphericWeight/`
	- `1-GenSimpleMC.py`: Generates a simple Monte Carlo sample.
	- `2-WriteGeneratorConfig.py`: Serializes LeptonInjector configuration and writes event tables.
	- `3-AtmFlux.py`: Calculates atmospheric neutrino flux using nuSQuIDS.
	- `4-GetWeights.py`: Calculates event weights using atmospheric flux.
	- `AtmFlux_output.h5`: Output from atmospheric flux calculation.
	- `config.lic`, `muons.i3`, `output.h5`, `nuweights.txt`: Configuration, input, and output files.

## Usage

1. Generate muon events with `1-GenSimpleMC.py`.
2. Write generator configuration and event tables with `2-WriteGeneratorConfig.py`.
3. For atmospheric weights, run `AtmosphericWeight/3-AtmFlux.py` to produce flux tables, then `AtmosphericWeight/4-GetWeights.py` to calculate weights.
4. For basic weights, run `BasicWeight/3-GetEventWeights.py`.


## Notes

- Atmospheric weighting uses nuSQuIDS for realistic flux calculations.
- Basic weighting uses a simple power-law flux.
- Input/output files are shared between steps as needed.

## Open Questions
How are and should multiple S frames handled?
Is the weighting correct?

