Calamari
========

This python project has tools to analyze SQUID data from the Berkeley DAQ.
It consists of a lightweight, modular framework to loop through ROOT trees with waveforms, apply triggering, filtering, calculations of parameters, etc.
It can then write the results to a ROOT file or a python pickle file

Dependencies
------------

Code Structure
--------------

The Calamari/ directory contains the core code.
The general container and controller class is the FileLooper (in base.py).
This class takes in a list of ROOT file names, and its methods loop over the tree holding the data and apply operations.
This class also takes in a Trigger module, as well as an arbitrary number of additional modules.
The Trigger module works on the ROOT tree entry directly, finding sections of raw waveform that satisfies some trigger condition.
When a waveform is triggered, the raw data is converted from ROOT format to a python dictionary, and passed to the additional modules.
Each module then acts on this dictionary, calculating quantities from the waveform, and adding the results to the dictionary
After finishing the loop, all the resulting data can be saved to an output file.

Currently implemented are two different triggers (in triggers.py).
The SimpleTrigger has a rolling time window, and triggers an event if the voltage in a given time window is a certain number of standard deviations above the average.
This was used to analyze Raul's first "background" run, where a number of thermal pulses of unknown (particle?) origin were observed.
The HeaterTrigger triggers the bolometer waveform when a heater pulse is seen. This requires two streams of data in the file - every other entry in the ROOT tree should contain a bolometer waveform readout or the heater pulser readout.
The HeaterTrigger searches the heater pulser readout for a pulse above a threshold, and when found it reads out the bolometer waveform data from the same time period.

Beyond these triggers, two additional modules are implemented (in filters.py).
The Filter class applies a basic Butterworth filter to a provide waveform, and returns the filtered waveform.
The PulseParams class calculates basic pulse parameters, such as amplitude, decay time, leading edge time, etc.

Example Scripts
---------------
