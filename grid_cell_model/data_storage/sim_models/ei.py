'''
Simulation data model for the EI network: in general a network of two
populations, E and I cells; more specifically, a grid cell network.
'''
import numpy as np

from analysis import spikes


def extractSpikes(mon):
    '''
    Extract spikes from a spike monitor (a dict-like object), that contains the
    relevant fields.
    
    Return a tuple (senders, spikeTimes).
    '''
    e = mon['events']
    return (e['senders'], e['times'])



class MonitoredSpikes(spikes.PopulationSpikes):
    '''
    This class extracts data from a DataStorage object (a dictionary).
    '''

    def __init__(self, data, monName, NName):
        '''
        Return the senders and spike times from a monitor in the data
        dictionary

        Parameters
        ----------
        data : dict
            A data dictionary containing the monitor monName.
        monName : string
            The spike monitor dictionary name.
        NName : string
            Name of the network parameter that specifies the total number of
            neurons.
        '''
        N = data['net_attr'][NName]
        senders, times = extractSpikes(data[monName])
        spikes.PopulationSpikes.__init__(self, N, senders, times)

