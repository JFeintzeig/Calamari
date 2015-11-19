import glob
import numpy
import pylab
import Calamari

######################################################
# Define file names, name of tree in root file, etc. #
######################################################

file_names=glob.glob('/Users/jfeintzeig/CUORE/Data/squid/Squid_Ch3_50to70mK_000027*.root')
tree_name='data_tree'

#############################################
# Create FileLooper, add trigger and filter #
#############################################

FL=Calamari.FileLooper(file_names,tree_name)
trigger=Calamari.SimpleTrigger('SimpleTrigger',0.01,1e-4,1e4,5,0.5)
FL.add_trigger(trigger)
butter=Calamari.Filter('Butterworth',30,5,1e-4,'Waveform','FilteredWaveform',True)
FL.add_module(butter)

#######
# Run #
#######

FL.loop(500)
FL.finish()
FL.write_output('simple_trigger_filter')
