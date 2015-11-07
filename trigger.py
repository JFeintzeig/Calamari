import ROOT
import numpy
from base import *

class Trigger(Module):
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
        super(Trigger,self).__init__(name)

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

    @Module._execute
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

    @Module._finish
    def finish(self):
        pass
