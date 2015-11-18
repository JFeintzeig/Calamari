import time
import numpy
import ROOT
import root_numpy

class FileLooper(object):
    '''
    Class to control program flow for analyzing SQUID data
    '''
    def __init__(self,file_names,tree_name,logging=True):
        '''

        Input:
        -file_names: list of strings corresponding to path of each ROOT
            file to load
        -tree_name: string with name of tree in ROOT file
        -logging: bool, set to True for periodic useful output
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
        '''
        Function to add trigger module

        Input:
        -trigger: trigger class instance
        '''
        self.trigger=trigger
        if self.logging:
            print "Added trigger %s" % (trigger.name)

    def add_module(self,module):
        '''
        Function to add modules to process events

        Input:
        -module: module class instance
        '''
        self.modules[module.name]=module
        if self.logging:
            print "Added module %s" % (module.name)

    def open_file(self):
        '''
        Creates ROOT TChain, adds files to chain
        '''
        self.chain=ROOT.TChain(self.tree_name)
        for f in self.file_names:
            self.chain.Add(f)
        if self.logging:
            print "Opened files"

    def loop(self,n_events=0):
        '''
        This is where all the actual calculations/work happens. This function
        loops through each entry of the TChain. For each entry, it:
            -applies the trigger module on the ROOT TChain entry, which returns
            a dict of data for that event
            -applies each additional module to this outputed dict. each
            module should return a new dict, containing previous and
            additional event data
            -saves final event dict in a list of dicts

        Input:
        -n_events: int, number of events to process. if 0, it processes all
        events in files
        '''
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
                # if an event doesn't trigger, trigger should return False
                # similarly, if a module rejects an event, module shoul return False
                if event:
                    event=self.modules[name].execute(event)
            if event:
                self.events+=[event]
            self.i+=1
            if self.logging and self.i%100==0:
                print "Event %i / %i" % (self.i, n_events)
            if self.i==n_events:
                break
        end_time=time.time()

    def write_output(self,outfile_name,ROOT=True,pickle=False):
        '''
        Converts outputed event data into re-usable data format,
        either ROOT or a pickled numpy record array

        Input:
        -outfile_name: string, name of outfile to write to (without extension)
        -ROOT: bool, True to write a ROOT file
        -pickle: bool, True to write a pickle file
        '''
        # build record array of outputed event dicts
        # assume all dict's are the same structure
        dt=[]
        for item in self.events[0]:
            if type(self.events[0][item])==numpy.ndarray:
                dt+=[(item,type(self.events[0][item][0]),len(self.events[0][item]))]
            else:
                dt+=[(item,type(self.events[0][item]))]

        dt=numpy.dtype(dt)
        values=[tuple(each.values()) for each in self.events]
        out=numpy.zeros((len(self.events),),dtype=dt)
        out[:]=values

        # convert record array to root tree, write to file
        if ROOT:
            if self.logging:
                print "Creating file %s.root" % (outfile_name)
            root_numpy.array2root(out,'%s.root' % (outfile_name))
        # write record array to pickle file
        if pickle:
            if self.logging:
                print "Creating file %s.pickle" % (outfile_name)
            f=open('%s.pickle' % (outfile_name),'w')
            pickle.dump(out,f)
            f.close()

    def finish(self):
        '''
        Calls finish methods of trigger and each modules
        '''
        self.trigger.finish()
        for name in self.modules:
            self.modules[name].finish()
        if self.logging:
            print "Finished"

class Module(object):
    '''
    Base class for Modules to plug into FileLooper

    All Modules can have names (strings), event counters, and timers

    For a Trigger Module, the execute() method should take in a TChain and an
    integer corresponding to the current position in the TChain to consider

    For a standard Module, the execute() method should take in a dict that contains
    all data for that given (already triggered) event
    '''
    def __init__(self,name):
        '''
        Input:
        -name: string name of module for nice printing later
        '''
        self.counter=0 # event counter
        self.name=name # string name
        self.run_time=0 # for timer

    @staticmethod
    def _execute(exec_fn):
        '''
        Static method to wrap another function in a timer and a counter. This can
        be added as a decorator to execute() methods of specific modules to easily
        add this functionality.
        '''
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
        '''
        Static method to wrap another function with a pretty printing display of
        timer results.  This can be added as a decorator to finish() methods of
        specific models to easily add this functionality.
        '''
        def wrapper(self,*args,**kwargs):
            finish_fn(self,*args,**kwargs)
            print type(self).__name__, "%s: %i events, %3.2f seconds" %\
                (self.name, self.counter, self.run_time)
        return wrapper

    def execute(self):
        '''
        This function is called on each event.  This is where the business happens.
        '''
        pass

    def finish(self):
        '''
        After processing all events, this function is called once.
        '''
        pass
