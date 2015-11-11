import sys
import time
import ROOT
import numpy
import pylab
from scipy.signal import butter, lfilter, freqz

class FileLooper(object):
    def __init__(self,file_names,tree_name,logging=True):
        '''
        file_names: list of strings corresponding to path to each file to load
        tree_name: string with name of tree in root file
        '''
        self.file_names=file_names
        self.tree_name=tree_name
        self.modules={}
        self.trigger=None
        self.events=[]
        self.logging=logging
        if self.logging:
            print "Created FileLooper"

    def add_trigger(self,trigger):
        self.trigger=trigger
        if self.logging:
            print "Added trigger %s" % (trigger.name)

    def add_module(self,module):
        self.modules[module.name]=module
        if self.logging:
            print "Added module %s" % (module.name)

    def open_file(self):
        self.chain=ROOT.TChain(self.tree_name)
        for f in self.file_names:
            self.chain.Add(f)
        if self.logging:
            print "Opened files"

    def loop(self,n_events=0):
        start_time=time.time()
        self.open_file()
        if self.logging:
            print "Starting loop"
        self.i=0
        if n_events==0:
            n_events=self.chain.GetEntries()
        while self.i < n_events:
            # apply trigger to each entry in root tree
            # this builds triggered waveforms, puts them in dict structure
            # other modules work with event data in dict form
            event=self.trigger.execute(self.chain,self.i)
            for name in self.modules:
                if event:
                    event=self.modules[name].execute(event)
            self.events+=[event]
            self.i+=1
            if self.logging and self.i%100==0:
                print "Event %i / %i" % (self.i, n_events)
            if self.i==n_events:
                break
        end_time=time.time()

    def write_output(self,outfile_name):
        #build record array of output from each module
        #then turn into rootfile
        pass

    def finish(self):
        self.trigger.finish()
        for name in self.modules:
            self.modules[name].finish()
        if self.logging:
            print "Finished"

class Module(object):
    '''
    base class for modules to plug into FileLooper
    '''
    def __init__(self,name):
        self.counter=0
        self.name=name
        self.run_time=0

    @staticmethod
    def _execute(exec_fn):
        def wrapper(self,event,*args,**kwargs):
            start_time=time.time()
            out=exec_fn(self,event,*args,**kwargs)
            self.counter+=1
            end_time=time.time()
            dt=end_time-start_time
            self.run_time+=dt
            return out
        return wrapper

    @staticmethod
    def _finish(finish_fn):
        def wrapper(self,*args,**kwargs):
            finish_fn(self,*args,**kwargs)
            print type(self).__name__, "%s: %i events, %3.2f seconds" %\
                (self.name, self.counter, self.run_time)
        return wrapper

    def execute(self):
        pass

    def finish(self):
        pass
