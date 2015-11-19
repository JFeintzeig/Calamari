import ROOT
import numpy
from base import *

class SimpleTrigger(Module):
    '''
    Module to loop through waveforms in chunks, trigger on pieces of waveforms
    that are a certain number of standard deviations above the mean baseline,
    and then return a larger chunk of the waveform around the triggered section
    as an "event". This is meant for "background" data, ie. data where every entry
    in the TChain is bolometer data, not interspersed with heater waveforms
    '''
    def __init__(self,name,dt_trigger,dt_sample,n_samples,n_sigma,readout_length):
        '''
        Inputs:
        -name: string name of module
        -dt_trigger: float, time interval to average together data, compare to
        previous baseline, and decide whether to trigger or not
        -dt_sample: float, time interval of one sample in the waveform
        -n_samples: int, number of samples in waveform in TChain
        -n_sigma: float, trigger threshold in number of standard deviations above
        baseline
        -readout_length: float, time in seconds corresponding to how much of the
        waveform around the trigger point to read out
        '''
        self.dt_trigger=dt_trigger
        self.dt_sample=dt_sample
        self.n_samples=int(n_samples)
        self.n_sigma=n_sigma
        # daq_buffer contains daq data from the previous two entries
        self.daq_buffer=numpy.zeros(2*n_samples)
        self.readout_length=readout_length
        self.stride=self.dt_trigger/self.dt_sample
        super(SimpleTrigger,self).__init__(name)

    #TODO: could probably cache means of strides to calc full mean more efficiently
    def trigger(self):
        '''
        Apply trigger condition

        Returns:
        -bool, True if trigger was triggered
        '''
        # mean and std calculated using the waveform buffer up to the current chunk
        mean=numpy.mean(self.daq_buffer[:-self.stride])
        std=numpy.std(self.daq_buffer[:-self.stride])
        return numpy.mean(self.daq_buffer[-self.stride:])>mean+self.n_sigma*std

    def build_event(self,chain,i,j):
        '''
        Extract a sub-piece of the waveform around the triggered chunk
        '''
        # put extracted waveform here
        waveform=numpy.zeros(self.readout_length/self.dt_sample)
        # indices of TChain waveform to start/end extraction
        start=int(j*self.stride-self.readout_length/2/self.dt_sample)
        end=int(j*self.stride+self.readout_length/2/self.dt_sample)
        if start>0 and end<=self.n_samples:
            waveform=numpy.array(chain.Waveform[start:end])
        elif start<0:
            # if readout start falls in previous waveform, go back 1 entry in chain
            waveform[-1*start:]=chain.Waveform[:end]
            chain.GetEntry(i-1)
            waveform[:-1*start]=chain.Waveform[self.n_samples+start:]
        elif end>=n_samples:
            # if readout end falls in next waveform, go forward 1 entry in chain
            waveform[:self.n_samples-start]=chain.Waveform[start:]
            chain.GetEntry(i+1)
            waveform[self.n_samples-start:]=chain.Waveform[:end-self.n_samples]

        # skip ahead in waveform to after readout window before once again looking
        # for chunks to trigger
        #TODO: eek, if readout_length, n_samples, and stride don't all match up,
        # this will cause stride to be shifted, may not be divisible into n_samples
        # evenly, and inner loop below (j) could try to read beyond end of waveform
        j_out=j+int(self.readout_length/2/self.dt_sample/self.stride)
        chain.GetEntry(i)
        return waveform, j_out

    def execute(self,chain,i):
        '''
        Do this on each entry in TChain.  Increments through waveform in chunks of
        length self.stride, comparing the average voltage in each chunk to the
        avg/standard deviation of the data leading up to the chunk, and decides to
        trigger and readout the given chunk and surrounding window, or not.

        Input:
        -chain: ROOT TChain holding SQUID data
        -i: int, current position in TChain
        '''
        chain.GetEntry(i)
        # we look at chunks in waveform of length self.stride
        # j is the index to loop through chunks of waveform
        j=0
        events=[]
        while j < self.n_samples/self.stride:
            # shift daq_buffer left, fill end with current stride
            self.daq_buffer[:-self.stride]=self.daq_buffer[self.stride:]
            start=int(self.stride*j)
            end=int(self.stride*(j+1))
            self.daq_buffer[-self.stride:]=chain.Waveform[start:end]
            # if current stride satisfies trigger condition, read out waveform
            # and form event
            if self.trigger():
                wf,j_out=self.build_event(chain,i,j)
                events+=[{'Waveform':wf}]
                j=j_out
                #TODO: add other event data (times, etc.) to dict!
                print "got event!"
            else:
                j+=1
        if len(events)==0:
            return False
        elif len(events)==1:
            return events[0]
        elif len(events)>1:
            #TODO: need to figure out how to put multiple events back into
            # control flow for FileLooper.loop() to deal with...
            return events[0]

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
