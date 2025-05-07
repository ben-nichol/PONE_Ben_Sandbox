# SALT - Radioactive Potassium Noise Module

This module generates photons from radioactive <sup>40</sup>K background events over a DOM in the Cascadia Basin. Given an I3 file with a tree of physics photons, this noise generator will create a new tree of photons from <sup>40</sup>K decay events within some time window around the physics hits.

## Module Parameters

`SEED`:            Seed for the random number generator.

`salinity`:        Salinity of seawater. Optional if you want to specify a particular salinity other than the known salinity of the Cascadia Basin.

`charFile`:        Path to the .pkl file holding the <sup>40</sup>K background characterization use by the generator.

`inputTreeName`:   Name of the tree in the input I3 file containing physics photons.

`outputTreeName`:  Name of the tree that will be added to the I3 file containing generated <sup>40</sup>K noise information.

`startPadding`:    Time padding to add before the physics photons for noise generation.

`endPadding`:      Time padding to add after the physics photons for noise generation.

`skipSingles`:     Choose whether to ignore <sup>40</sup>K events that only hit the DOM with 1 photon.

`pmtHitsOnly`:     Choose whether to ignore all photons that didn't hit a PMT.

`skipOneFold`:     Choose whether to ignore any events that only deposit photons into one PMT. If this is set to True, pmtHitsOnly must also be set to True.

More information on each parameter is available in the module code.

## ToDo
- [ ] Add quantum efficiency cuts to generated photons
- [ ] Improve model for nearby events that form Cherenkov rings on the DOM
- [ ] Check and verify module output after incorporation into pone_offline

## Files

`Salt.py`

The <sup>40</sup>K noise generator module written in Python. Given an input I3 file with existing htis, this module will add I3 frames to the file filled with <sup>40</sup>K hits. It will generate as many frames as there are in the referenced tree of the passed I3 file.

----
`NoiseCharacterization.pkl`

<sup>40</sup>K noise signal characterization stored as a pickle file. This object holds distributions and functions for generating the noise signal.

----
`NoiseGenHelpers.py`

Python file that holds helpful functions and classes for the `Salt.py` module.

----
## `Examples` Directory
`IcetrayTest.py`

Sample Icetray script used to test the execution of the noise module.

---
`makeK40Hits.sh`

Bash script used to run `IcetrayTest.py` in batch mode.

---
`submit_to_condor.sh`

Bash script used to submit to condor.

---

If you have any questions contact me on Slack or email - Jakub Stacho - jakubs@sfu.ca
