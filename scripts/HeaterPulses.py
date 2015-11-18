import glob
import ROOT
import numpy
import root_numpy
import pylab
import Calamari

sec = 1.
msec = 1e-3*sec
window_time=30*msec
#trigger_voltage_thresh=-0.015
trigger_voltage_thresh=-0.025

#data_dir='/Users/jfeintzeig/CUORE/Data/squid/Run_000077/'
#out_dir='/Users/jfeintzeig/CUORE/Scripts/2015/SQUIDDataAnalysis/plots/Run_000077/'
data_dir='/Users/jfeintzeig/CUORE/Data/squid/Run_000103/'
out_dir='/Users/jfeintzeig/CUORE/Scripts/2015/SQUIDDataAnalysis/plots/Run_000103/'
file_names=glob.glob(data_dir+'*.root')
file_names=glob.glob(data_dir+'*.root')
tree_name='data_tree'

# create looper, add trigger and module
FL=Calamari.FileLooper(file_names,tree_name)
trigger=Calamari.HeaterTrigger('HeaterTrigger',window_time,trigger_voltage_thresh)
FL.add_trigger(trigger)
pulse_params=Calamari.PulseParams('PulseParams','Waveform')
FL.add_module(pulse_params)

# run!
FL.loop()
FL.finish()

#plot stuff
import pandas
events=[]
for item in FL.events:
    if item:
        events+=[item]

data=pandas.DataFrame(events)

sys.exit()

# for run 103
plot_vars=['HeaterEnergy','HeaterWidth','Time','PulseParams_DecayTime',\
    'HeaterAmplitude','PulseParams_PulseHeight','PulseParams_Baseline']

t0=1.4472008e9+265
from scipy import stats
width=data.HeaterWidth==0.00002
high=data.HeaterEnergy>0.4e-7
low=data.HeaterEnergy<0.4e-7
x=data['HeaterEnergy'][width]
y=data['PulseParams_PulseHeight'][width]
m1,b1,r,p,std=stats.linregress(x[high],y[high])
m2,b2,r,p,std=stats.linregress(x[low],y[low])
line=lambda x,m,b: m*x+b

pylab.plot(x[high],line(x[high],m1,b1),color='red',linewidth=1.3)
#pylab.plot(x[low],line(x[low],m2,b2),color='red',linewidth=1.3)
pylab.scatter(x,y,color='black')
pylab.xlabel(r'Heater "Energy" (V$^2$*s)')
pylab.ylabel('Pulse Height (V)')
pylab.xlim(0,1.7e-7)
pylab.ylim(0,0.025)
pylab.grid()
pylab.savefig(out_dir+'height_vs_energy.png')
pylab.show()

data.DT=data.PulseParams_LeadingEdgeTime-data.HeaterLeadingEdgeTime
good=data.DT>-0.2
pylab.hist(1000*data.DT[good],
    300,histtype='step',linewidth=2)
pylab.xlabel(r'$\Delta t$ (ms)',size='14')
pylab.ylabel('Counts')
pylab.grid()
pylab.xlim(0,1)
pylab.savefig(out_dir+'time_resolution.png')
pylab.show()

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
    #pylab.scatter(data[var],\
    #    data['PulseParams_PulseHeight'])
    pylab.xlabel(var)
    pylab.ylabel('PulseParams_PulseHeight')
    pylab.grid()
    pylab.show()

def timeme(fn,n,*args):
    i=0
    st=time.time()
    while i<n:
        fn(*args)
        i+=1
    et=time.time()
    print et-st

#TODO

#plots:
#-time diff between heater and bolo rising edge
#-waveforms of a few pulses with varying heights
#-waveforms of double-pulse events

#TODO: raul's RT curve plot

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

time_cut=data['Time']<t0
pylab.scatter(data['HeaterWidth'][time_cut],\
    data['PulseParams_PulseHeight'][time_cut],color='black')
pylab.xlim(0.0001,0.0008)
pylab.ylim(0,0.035)
pylab.xlabel('Heater Width (s)')
pylab.ylabel('Pulse Height (V)')
pylab.grid()
pylab.savefig(out_dir+'height_vs_heaterwidth.png')
pylab.show()

pylab.scatter(data['HeaterEnergy'][time_cut],\
    data['PulseParams_PulseHeight'][time_cut],color='black')
pylab.xlabel(r'Heater "Energy" (V$^2$*s)')
pylab.ylabel('Pulse Height (V)')
pylab.xlim(0.5e-7,4e-7)
pylab.ylim(0,0.035)
pylab.grid()
pylab.savefig(out_dir+'height_vs_heaterenergy.png')
pylab.show()

t=numpy.linspace(0,len(data.Waveform.loc[0])*1e-5,len(data.Waveform.loc[0]))
pylab.plot(t,data.Waveform.loc[0],linewidth=2,color='black')
pylab.grid()
pylab.xlabel('Time (s)')
pylab.ylabel('Voltage (V)')
pylab.savefig(out_dir+'example_pulse.png')
pylab.show()

t=numpy.linspace(0,len(data.Waveform.loc[0])*1e-5,len(data.Waveform.loc[0]))
for event in events[:15]:
    waveform=event['Waveform']
    ind=numpy.argmax(waveform)
    baseline=numpy.mean(waveform[0:int(0.75*ind)])
    std=numpy.std(waveform[0:int(0.75*ind)])
    leading_edge=min(numpy.argwhere(waveform>baseline+5*std))[0]*dt_sample+\
        event['Time']
    pylab.plot(event['Time']+t,waveform)
    pylab.axvline(leading_edge)
    pylab.show()

# small pulse
dts=1e-5
pylab.clf()
pylab.subplot(211)
f.data_tree.GetEntry(28)
t=numpy.linspace(0,1,1e5)
pylab.title('Bolometer')
pylab.plot(t,f.data_tree.Waveform,linewidth=2,color='black')
pylab.ylabel('Voltage (V)')
pylab.xlim(500*dts,5000*dts)
pylab.ylim(-0.035,-0.032)
pylab.subplot(212)
f.data_tree.GetEntry(29)
pylab.title('Heater')
pylab.plot(t,f.data_tree.Waveform,linewidth=2,color='black')
pylab.ylabel('Voltage (V)')
pylab.xlabel('Time (s)')
pylab.xlim(500*dts,5000*dts)
pylab.ylim(-0.027,0.002)
pylab.subplot(212)
pylab.savefig(out_dir+'small_pulse.png')
pylab.show()

'''
chain=ROOT.TChain("data_tree")
chain.Add(data_dir+"*.root")

events=[]
i=0
while i < 3:
    e=trigger.execute(chain,i)
    if e:
        pulse_params.execute(e)
        pylab.plot(e['Waveform'])
        pylab.title('%i %i' % (i,numpy.argmax(e['Waveform'])))
        pylab.show()
        events+=[e]
    i+=1
'''
