import numpy
import pylab
from scipy.signal import butter, lfilter, freqz
from base import *

class Filter(Module):
    '''
    apply basic Butterworth filter to waveform
    '''
    #TODO: have some dt_sample, i, num lieing around, change this...
    def __init__(self,name,cutoff,order,input_name, output_name):
        self.cutoff=cutoff
        self.order=order
        self.input_name=input_name
        self.output_name=output_name
        super(Filter,self).__init__(name)

    def butter_lowpass(self,fs):
        nyq = 0.5 * fs
        high = self.cutoff / nyq
        b, a = butter(self.order, high, btype='low')
        return b, a

    def butter_lowpass_filter(self,data,fs):
        b, a = butter_lowpass(fs,)
        y = lfilter(b, a, data)
        return y

    def plot_filtered(self):
        pylab.plot(t,wf,linewidth=2,color='black',label='Unfiltered Pulse')
        pylab.plot(t,filtered,linewidth=2,color='red',label='Filtered Pulse')
        pylab.title('Event number (Time in run in seconds): %i' % (i))
        pylab.grid()
        pylab.xlabel('Time (ms)')
        pylab.savefig(out_dir+'pulse_%i_%i.png' % (i,num))
        pylab.clf()

    @Module._execute
    def execute(self,event):
        wf=event[self.input_name]
        t=1000*numpy.linspace(0,len(wf)*dt_sample-dt_sample,len(wf))
        filtered=butter_lowpass_filter(wf,1/dt_sample)
        self.plot_filtered()
        event[self.output_name]=filtered
        return event

    @Module._finish
    def finish(self):
        pass

class PulseParams(Module):
    '''
    calculate pulse height, decay constant, etc.
    '''
    def __init__(self,name,waveform_name):
        self.waveform_name=waveform_name
        super(PulseParams,self).__init__(name)

    @Module._execute
    def execute(self,event):
        dt_sample=event['dt_sample']
        waveform=event[self.waveform_name]

        # calculate height
        ind=numpy.argmax(waveform)
        baseline=numpy.mean(waveform[0:int(0.75*ind)])
        std=numpy.std(waveform[0:int(0.75*ind)])
        height=waveform[ind]-baseline

        # calculate decay time
        thresh=baseline+height*(1/numpy.e)
        try:
            decay_time=(max(numpy.argwhere(waveform>thresh))[0]-ind)*dt_sample
        #TODO: had problems with pulses that are just flat lines, no data...
        # why does this happen??
        # fix this/treat this better...
        except ValueError:
            return False

        #calculate leading edge
        try:
            leading_edge=min(numpy.argwhere(waveform>baseline+5*std))[0]*dt_sample+\
                event['Time']
        #TODO: if no samples 5 sigma above baseline...need to treat this better
        except ValueError:
            leading_edge=-1

        # write output
        event[self.name+'_LeadingEdgeTime']=leading_edge
        event[self.name+'_Baseline']=baseline
        event[self.name+'_PulseHeight']=height
        event[self.name+'_DecayTime']=decay_time
        return event

    @Module._finish
    def finish(self):
        pass

