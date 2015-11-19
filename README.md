Calamari
========

This python project has tools to analyze SQUID data from the Berkeley DAQ.
It consists of a lightweight, modular framework to loop through ROOT trees with waveforms, apply triggering, filtering, calculations of parameters, etc.
It can then write the results to a ROOT file or a python pickle file

Installation
------------

No installation is required.  Just clone the repository, add the directory to your $PYTHONPATH, and you're good to go:

>In [1]: import Calamari

>In [2]: Calamari.
>Calamari.FileLooper(        Calamari.__getattribute__(  Calamari.__subclasshook__(
>Calamari.Filter(            Calamari.__hash__(          Calamari.base
>Calamari.HeaterTrigger(     Calamari.__init__(          Calamari.butter(
>Calamari.Module(            Calamari.__name__           Calamari.filter
>Calamari.PulseParams(       Calamari.__new__(           Calamari.freqz(
>Calamari.ROOT               Calamari.__package__        Calamari.lfilter(
>Calamari.SimpleTrigger(     Calamari.__path__           Calamari.numpy
>Calamari.__class__(         Calamari.__reduce__(        Calamari.pylab
>Calamari.__delattr__(       Calamari.__reduce_ex__(     Calamari.root\_numpy
>Calamari.__dict__           Calamari.__repr__(          Calamari.time
>Calamari.__doc__            Calamari.__setattr__(       Calamari.trigger
>Calamari.__file__           Calamari.__sizeof__(        
>Calamari.__format__(        Calamari.__str__(           


Requirements and Dependencies
-----------------------------

This project was created using:

>Python 2.7.9

>ROOT 5.34/25 w/pyROOT pybindings

The following python libraries are required for processing data:

>numpy==1.9.2

>root-numpy==4.4.0

>scipy==0.15.1

The following libraries are helpful for plotting, and some of the example scripts may make use of them:

>matplotlib==1.4.2

>pandas==0.16.2

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

In the scripts/ directory, the HeaterPulses.py script shows how to process some basic heater pulse bolometer data.
This is the script I used to analyze the pulser amplitude scan data and plot the linearity and time resolution of the bolometer.

The TestSimpleTrigger.py script shows how to apply the SimpleTrigger and a Butterworth filter to "background" data (ie. data that is not triggered by heater pulses).
