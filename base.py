import sys
import ROOT
import numpy
import pylab
from scipy.signal import butter, lfilter, freqz

class FileLooper(object):
    def __init__(self,file_names,tree_name):
        '''
        file_names: list of strings corresponding to path to each file to load
        tree_name: string with name of tree in root file
        '''
        self.file_names=file_names
        self.tree_name=tree_name
        self.modules={}
        self.trigger=None

    def add_trigger(self,trigger):
        self.trigger=trigger

    def add_module(self,module):
        self.modules[module.name]=module

    def open_file(self):
        self.chain=ROOT.TChain(self.tree_name)
        for f in self.file_names:
            self.chain.Add(f)

    def loop(self,n_events=0):
        self.open_file()
        i=0
        if n_events==0:
            n_events=self.chain.GetEntries()
        while i < n_events:
            # apply trigger to each entry in root tree
            # this builds triggered waveforms, puts them in dict structure
            # other modules work with event data in dict form
            events=self.trigger.execute(self.chain,i)
            for event in events:
                for name in self.modules:
                    self.modules[name].execute(event)
            i+=1
            if i==n_events:
                break

    def write_output(self,outfile_name):
        #build record array of output from each module
        #then turn into rootfile
        rec
        root_numpy.array2root(rec,outfile_name)

    def finish(self):
        for name in self.modules:
            self.modules[name].finish()

# modules store the info they calculate?
# how to put info back into tree??
class Module(object):
    '''
    base class for modules to plug into FileLooper
    '''
    def __init__(self,name):
        self.counter=0
        self.name=name

    @staticmethod
    def _execute(exec_fn):
        def wrapper(self,event):
            exec_fn(self,event)
            self.counter+=1
        return wrapper

    @staticmethod
    def _finish(finish_fn):
        def wrapper(self):
            finish_fn(self)
            print type(self).__name__, "%s: %i events" % (self.name, self.counter)
        return wrapper

    def execute(self):
        pass

    def finish(self):
        pass

class Average(Module):
    '''
    simple module to calculate the average baseline for an event
    '''
    def __init__(self,name,input_name):
        self.avg=[]
        self.input_name=input_name
        super(Average,self).__init__(name)

    @Module._execute
    def execute(self,event):
        avg=numpy.mean(event[self.input_name])
        print avg
        self.avg+=[avg]

    @Module._finish
    def finish(self):
        self.avg=numpy.array(self.avg)
