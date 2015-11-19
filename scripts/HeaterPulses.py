import glob
import numpy
import pylab
import Calamari

###########################################
# Setup various constants and directories #
###########################################

SEC = 1.
MSEC = 1e-3*SEC

# when we find heater pulses, read out a 30 ms window around the pulse
window_time=30*MSEC
# only trigger on heater pulses with voltages beyond this number
trigger_voltage_thresh=-0.025 # Volts

# this directory should contain only ROOT files from this specific run
# when I ran it I only used the first ~120 files
data_dir='/Users/jfeintzeig/CUORE/Data/squid/Run_000103/test/'
# directory to save plots in
out_dir='/Users/jfeintzeig/CUORE/Scripts/2015/SQUIDDataAnalysis/plots/Run_000103/test/'
file_names=glob.glob(data_dir+'*.root')
file_names=glob.glob(data_dir+'*.root')
# name of tree in root file where data is stored
tree_name='data_tree'

########################################################################
# Setup FileLooper class, triggers, modules, and put them all together #
########################################################################

FL=Calamari.FileLooper(file_names,tree_name)

trigger=Calamari.HeaterTrigger('HeaterTrigger',window_time,trigger_voltage_thresh)
FL.add_trigger(trigger)

pulse_params=Calamari.PulseParams('PulseParams','Waveform')
FL.add_module(pulse_params)

#####################
# Process the data! #
#####################

FL.loop()
FL.finish()
FL.write_output('heater_pulses')

##################################
# Extra stuff - make plots, etc. #
##################################

# Much of this stuff is specific to Run 103, the first run with heater pulse
# height scans that I looked at.  I include it below just as an example of how
# to work with this data and make plots in python.  This can all easily be done
# in ROOT instead as well

plot_heater_data=False
if plot_heater_data:

    import pandas
    data=pandas.DataFrame(FL.events)

    # select events with a certain heater width
    # fit pulse height vs. heater energy with a linear function
    t0=1.4472008e9+265
    from scipy import stats
    width=data.HeaterWidth==0.00002
    high=data.HeaterEnergy>0.4e-7 # only fit high-energy data
    x=data['HeaterEnergy'][width]
    y=data['PulseParams_PulseHeight'][width]
    m1,b1,r,p,std=stats.linregress(x[high],y[high])
    line=lambda x,m,b: m*x+b

    # plot pulse height vs. heater energy, along with linear fit
    pylab.plot(x[high],line(x[high],m1,b1),color='red',linewidth=1.3)
    pylab.scatter(x,y,color='black')
    pylab.xlabel(r'Heater "Energy" (V$^2$*s)')
    pylab.ylabel('Pulse Height (V)')
    pylab.xlim(0,1.7e-7)
    pylab.ylim(0,0.025)
    pylab.grid()
    pylab.savefig(out_dir+'height_vs_energy.png')
    pylab.show()

    # plot time resolution between heater leading edge and bolometer
    # pulse leading edge
    data.DT=data.PulseParams_LeadingEdgeTime-data.HeaterLeadingEdgeTime
    # remove large negative outliers, which are likely due to the module
    # failing to find a good leading edge time on the bolometer pulse
    good=data.DT>-0.2
    pylab.hist(1000*data.DT[good],
        300,histtype='step',linewidth=2)
    pylab.xlabel(r'$\Delta t$ (ms)',size='14')
    pylab.ylabel('Counts')
    pylab.grid()
    pylab.xlim(0,1)
    pylab.savefig(out_dir+'time_resolution.png')
    pylab.show()

    # make histograms and scatter plots of these variables
    plot_vars=['HeaterEnergy','HeaterWidth','Time','PulseParams_DecayTime',\
        'HeaterAmplitude','PulseParams_PulseHeight','PulseParams_Baseline']
    for var in plot_vars:
        pylab.hist(data[var],100,histtype='step',color='black',linewidth=2)
        pylab.xlabel(var)
        pylab.grid()
        pylab.show()
        pylab.scatter(data['Time'],data[var])
        pylab.xlabel('Time')
        pylab.ylabel(var)
        pylab.grid()
        pylab.show()
        pylab.scatter(data[var][data['Time']<t0],\
            data['PulseParams_PulseHeight'][data['Time']<t0])
        pylab.xlabel(var)
        pylab.ylabel('PulseParams_PulseHeight')
        pylab.grid()
        pylab.show()

    # pulse height histogram and decay time histogram for a given heater energy
    heater_energy=(data['HeaterEnergy']>1.25e-7)&(data['HeaterEnergy']<1.4e-7)
    data['PulseParams_PulseHeight'][heater_energy].hist(bins=40,
        color='black',histtype='step',linewidth=2,grid=True)
    pylab.xlabel('Pulse Height (V)')
    pylab.xlim(0.005,0.030)
    pylab.savefig(out_dir+'pulse_height.png')
    pylab.show()

    data['PulseParams_DecayTime'][heater_energy].hist(bins=40,
        color='black',histtype='step',linewidth=2,grid=True)
    pylab.xlabel('Decay Time (s)')
    pylab.xlim(0,0.002)
    pylab.savefig(out_dir+'decay_time.png')
    pylab.show()
