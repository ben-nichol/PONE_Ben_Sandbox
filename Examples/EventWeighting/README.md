# EventWeighting Example Scripts

This directory contains example scripts for calculating the weights for MC. There are two subdirectories which show two different approaches.

## Input
In both subdirectories an example set of input can be produced with:
- `1-GenSimpleMC.py`: Generates 10 muons between 100-1000GeV.
- `2-WriteGeneratorConfig.py`: Reads one or more .i3 files to write the LeptonInjector config.lic and data tables. You can specify any number of input files, followed by the output file name (used for both .h5 and .lic outputs). Example usage:
		python3 2-WriteGeneratorConfig.py input1.i3 input2.i3 ... outputName

## Weighting
Each subdirectory gives an option for generating weights:

- `BasicWeight/`
	- `3-GetEventWeights.py`: Calculates event weights using a basic power-law flux. Supports multiple LIC and HDF5 input files; all generators are combined and each event is weighted against all generators. Example usage:
	
	    python3 3-GetEventWeights.py --lic config1.lic config2.lic --input data1.h5 data2.h5
	
	  This will process all specified LIC and HDF5 files and print event weights for each event.

- `AtmosphericWeight/`
	- `3-AtmFlux.py`: Calculates atmospheric neutrino flux using nuSQuIDS.
	- `4-GetWeights.py`: Calculates event weights using atmospheric flux. Supports multiple LIC and HDF5 input files; all generators are combined and each event is weighted against all generators. Example usage:
	
	    python3 4-GetWeights.py --lic config1.lic config2.lic --events data1.h5 data2.h5 --nusquids_flux AtmFlux_output.h5
	
	  This will process all specified LIC and HDF5 files and write event weights to the output file.


## Open Questions
Is the weighting correct?
