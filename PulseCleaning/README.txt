The PulseCleaning folder contains methods used to clean up the pulse series to remove noise and stray light pulses that mess up the reconstructions.

Files:
	SignificantHitPulseCleaning.py - Removes pulses that are not within a window, if the DOM doesn't see enought light then all pulses are kept. 

	RecoPulseGenerator.py - Depricated method for converting MCPhotons to RecoPulses.

	TimeShift.py - A nearly depricated method used to shift the times in the event to a more reasonable time.
