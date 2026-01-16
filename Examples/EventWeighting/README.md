
# EventWeighting Example Scripts

This directory contains example scripts for calculating the weights for MC. There are two subdirectories which show two different approaches

## Structure
- `BasicWeight/`
	- `1-GenSimpleMC.py`: Generates 10 muons between 100-1000GeV.
	- `2-WriteGeneratorConfig.py`: Reads the .i3 file to write the LeptonInjector config.lic and data tables.
	- `3-GetEventWeights.py`: Calculates event weights using a basic power-law flux.

- `AtmosphericWeight/`
	- `1-GenSimpleMC.py`: Generates 10 muons between 100-1000GeV.
	- `2-WriteGeneratorConfig.py`: Reads the .i3 file to write the LeptonInjector config.lic and data tables.
	- `3-AtmFlux.py`: Calculates atmospheric neutrino flux using nuSQuIDS.
	- `4-GetWeights.py`: Calculates event weights using atmospheric flux.

## Usage

1. Generate muon events with `1-GenSimpleMC.py`.
2. Write generator configuration and event tables with `2-WriteGeneratorConfig.py`.
3. For atmospheric weights, run `AtmosphericWeight/3-AtmFlux.py` to produce flux tables, then `AtmosphericWeight/4-GetWeights.py` to calculate weights.
4. For basic weights, run `BasicWeight/3-GetEventWeights.py`.


## Open Questions
How are and should multiple S frames be handled?
Is the weighting correct?
