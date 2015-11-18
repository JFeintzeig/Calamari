import ROOT
import numpy
from base import *

class SimpleTrigger(Module):
    '''
    module to loop through waveforms in chunks, trigger on pieces of waveforms
    that are a certain number of standard deviations above the mean baseline,
    and then return a larger chunk of the waveform around the triggered section
    as an "event"

    while most other Modules should take in an event dictionary, this module
    instead takes in the TChain and the current entry number to look at...
    '''
    def __init__(self,name,dt_trigger,dt_sample,n_samples,n_sigma):
        self.dt_trigger=dt_trigger
        self.dt_sample=dt_sample
        self.stride=self.dt_trigger/self.dt_sample
        self.n_samples=n_samples
        self.n_sigma=n_sigma
        self.daq_buffer=numpy.zeros(2*n_samples)
        super(SimpleTrigger,self).__init__(name)

    #TODO: i just copy/pasted the following 31 lines, need to fix up for module,
    # integrate with standard module stuff...
    #TODO: could probably cache means of strides to calc full mean more efficiently
    def trigger(self):
        mean=numpy.mean(self.daq_buffer[:-self.stride])
        std=numpy.std(self.daq_buffer[:-self.stride])
        return numpy.mean(self.daq_buffer[-self.stride:])>mean+self.n_sigma*std

    def build_event(self,chain,i,j):
        #readout 0.5 sec interval, 0.25 sec before trigger point, 0.25 sec after
        waveform_length=0.5
        waveform=numpy.zeros(0.5/self.dt_sample)
        start=int(j*self.stride-0.25/self.dt_sample)
        end=int(j*self.stride+0.25/self.dt_sample)
        if start>0 and end<=self.n_samples:
            waveform=numpy.array(chain.Waveform[start:end])
        elif start<0:
            waveform[-1*start:]=chain.Waveform[:end]
            chain.GetEntry(i-1)
            waveform[:-1*start]=chain.Waveform[self.n_samples+start:]
        elif end>=n_samples:
            waveform[:self.n_samples-start]=chain.Waveform[start:]
            chain.GetEntry(i+1)
            waveform[self.n_samples-start:]=chain.Waveform[:end-self.n_samples]

        #TODO: fix this
        #eek, if waveform_length, n_samples, and stride don't all match up,
        # this will cause stride to be shifted, may not be divisible into n_samples
        # evenly, and inner loop below (j) could try to read beyond end of waveform
        j_out=j+int(waveform_length/2/self.dt_sample/self.stride)
        chain.GetEntry(i)
        return waveform, j_out

    def execute(self,chain,i):
        # here we do the inner loop over j's, calling build_event for each...
        chain.GetEntry(i)
        j=0
        events=[]
        while j < n_samples/stride:
            self.daq_buffer[:-self.stride]=self.daq_buffer[self.stride:]
            start=int(self.stride*j)
            end=int(self.stride*(j+1))
            self.daq_buffer[-self.stride:]=chain.Waveform[start:end]
            #trigger
            if self.trigger():
                #get waveform of event
                wf,j_out=self.build_event(chain,i,j)
                #TODO: change this to build event from all info in chain, turn into dict
                events+=[wf]
                j=j_out
                print "got event!"
            else:
                j+=1
        return events

    def finish(self):
        pass

class HeaterTrigger(Module):
    '''
    Trigger to find bolometer pulses from heater events
    For a given tree entry that contains data from Ch 0 (aka the bolometer),
    this expects the previous tree entry to contain data from Ch 1 (the heater pulser)
    taken at the exact same time.  It uses the heater data from the previous entry to
    search for pulses above a certain voltage threshold.  If it finds such a pulse,
    it finds the trigger time and returns the bolometer data from that time and a
    surrounding time window
    '''
    def __init__(self,name,window_time,min_voltage_trigger_threshold):
        '''
        Input:
        -name: string, name to call module for pretty printing
        -window_time: float, time in seconds to read-out around pulse. this will be
        the total time of the waveform extracted from each event
        -min_voltage_trigger_threshold: float, the heater voltage threshold below
        which we decide the data has a pulse and trigger. this is set-up to assume
        the heater pulses are negative (neg. spikes in voltage) while the bolometer
        pulses show a baseline increase
        '''
        self.window_time=window_time
        self.min_voltage_trigger_threshold=min_voltage_trigger_threshold
        # these will be filled from the first event
        self.window=0
        self.dt_sample=0
        self.n_samples=0
        super(HeaterTrigger,self).__init__(name)

    def calculate_heater_params(self,chain_entry,min_ind):
        '''
        Function to calculate parameters of the heater pulse. This could also
        be moved to a separate module in the future, since it doesn't really
        depend on the trigger.  But for now its convenient to have it here.

        Input:
        -chain_entry: entry in ROOT TChain corresponding to a piece of
        heater pulse data
        -min_ind: index of Waveform array corresponding to the minimum (ie. the
        peak of the negative heater pulse)

        Returns:
        -float, the amplitude of the heater pulse, in volts
        -float, the width of the heater pulse, in seconds
        -float, the energy proxy of the heater pulse, calculated as the
        voltage^2 * width
        -float, the leading edge time of the heater pulse, in seconds from the
        waveform start

        '''
        w=numpy.array(chain_entry.Waveform)
        # width is defined as the amount of time the square heater pulse is beyond
        # 90% of its peak value.  remember: the heater pulse is a negative spike in
        # voltage
        thresh=w[min_ind]*0.9
        pulse=w[w<thresh]
        heater_width=len(pulse)*chain_entry.dt_s # seconds

        heater_leading_edge_time=min(numpy.argwhere(w<thresh))[0]*self.dt_sample
        # baseline is mean value of waveform, excluding samples that have 90%
        # or more of the heater pulse peak voltage
        baseline=numpy.mean(w[w>thresh]) # mean value of waveform, excluding
        heater_amp=numpy.mean(pulse)-baseline
        # energy proxy, need resistance to convert to energy
        heater_energy=heater_amp**2*heater_width
        return heater_amp, heater_width, heater_energy, heater_leading_edge_time

    def read_out_waveform(self,chain,i,ind):
        '''
        Isolates a sub-section of the total waveform that is centered on the pulse

        Input:
        -chain: TChain containing heater and bolometer data
        -i: TChain index corresponding to current bolometer event to trigger
        -ind: index in waveform array corresponding to peak of heater pulse

        Returns:
        -waveform array containing bolomoter data at time triggered by heater pulse,
        +/- a time window given by total time window self.window_time
        -int corresponding to the index in the array where the readout window will
        start.  can be negative (gives position in waveform of previous bolo event)
        '''
        chain.GetEntry(i)
        start=ind-self.window/2
        end=ind+self.window/2
        waveform=numpy.zeros(end-start)
        if start<0:
            # if readout window spills into previous pulse, get that pulse from chain
            waveform[-1*start:]=chain.Waveform[:end]
            # chain alternates heater/bolometer entries, so go back 2 entries
            # to get last bolometer data
            chain.GetEntry(i-2)
            waveform[:-1*start]=chain.Waveform[self.n_samples+start:]
        elif end>self.n_samples:
            # if readout window spills into next pulse
            waveform[:self.n_samples-start]=chain.Waveform[start:]
            chain.GetEntry(i+2) # again, jump 2 spots
            waveform[self.n_samples-start:]=chain.Waveform[:end-self.n_samples]
        else:
            # readout window fully contained in current waveform
            waveform=numpy.array(chain.Waveform[start:end])
        return waveform, start

    @Module._execute
    def execute(self,chain,i):
        '''
        Execute this function on every entry in chain. It checks for a heater pulse
        beyond a threshold, and if found, reads out bolometer data around that time.

        Input:
        -chain: ROOT TChain containing data to process
        -i: int corresponding to current position in the TChain

        Output:
        -dictionary containing triggered bolometer waveform, heater waveform for same
        time window, time of the start of the window, and various heater pulse
        parameters
        '''
        event={}
        chain.GetEntry(i)

        # ch 1 is heater channel, i want bolometer events
        if chain.Ch==1:
            return False
        assert chain.Ch==0

        # fill info that should be the same for all events
        if self.n_samples==0:
            self.n_samples=len(chain.Waveform)
        if self.dt_sample==0:
            self.dt_sample=chain.dt_s
            self.window=int(self.window_time/self.dt_sample)

        # get event time
        t_s=chain.t_s
        t_mus=chain.t_mus
        event['Time']=t_s+t_mus*1e-6
        event['dt_sample']=chain.dt_s

        # get corresponding heater pulse
        chain.GetEntry(i+1)
        assert chain.t_s==t_s
        assert chain.t_mus==t_mus
        ind=numpy.argmin(chain.Waveform)

        # throw away pulses below heater trigger threshold
        if chain.Waveform[ind]>self.min_voltage_trigger_threshold:
            return False

        # find time,energy,etc of heater pulse
        heater_amp, heater_width, heater_energy, heater_leading_edge_time=\
            self.calculate_heater_params(chain,ind)
        event['HeaterAmplitude']=heater_amp
        event['HeaterWidth']=heater_width
        event['HeaterEnergy']=heater_energy
        event['HeaterLeadingEdgeTime']=heater_leading_edge_time+event['Time']

        # read out bolometer waveform for heater time +/- window/2
        event['Waveform'],start=self.read_out_waveform(chain,i,ind)
        event['HeaterWaveform'],start=self.read_out_waveform(chain,i+1,ind)

        event['Time']+=start*self.dt_sample

        # return event dictionary
        return event

    @Module._finish
    def finish(self):
        pass
