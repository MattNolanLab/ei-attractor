'''
Simulation data model for the EI network: in general a network of two
populations, E and I cells; more specifically, a grid cell network.
'''
import numpy as np

from ...analysis import spikes


def getNetParam(data, p):
    '''Extract a network parameter (p) from the data dictionary'''
    return data['net_attr'][p]

def getOption(data, o):
    return data['options'][o]

def extractSpikes(mon):
    '''
    Extract spikes from a spike monitor (a dict-like object), that contains the
    relevant fields.

    Return a tuple (senders, spikeTimes).
    '''
    e = mon['events']
    return (e['senders'], e['times'])


def sliceSignal(t, sig, tStart, tEnd):
    idx = np.logical_and(t >= tStart, t <= tEnd)
    return t[idx], sig[idx], idx


def extractStateVariable(mon, nIdx, varStr):
    '''Extract state variable from a monitor.

    Parameters
    ----------
    mon : list of dicts
        A list of (NEST) monitors, each monitoring one neuron.
    nIdx : int
        Neuron index
    varStr : str
        Name of the variable
    output
        A tuple (data, dt), for the signal
    '''
    n = mon[nIdx]
    return n['events'][varStr], n['interval']


def sumAllVariables(mon, nIdx, varList):
    '''
    Extract all variables from the list of monitors and sum them. The variables
    must implement the + operator.

    Parameters
    ----------
    mon : a list of dicts
        A list that contains dictionaries of monitors. The list should be
        compatible with the extractStateVariable function.
    nIdx : int
        Neuron index
    varList : list of strings
        Contains the list of variables that whould be extracted from the
        monitor and summed up.
    output
        A tuple (sum, dt) that contains the sum of all the variables 'sum' and
        the sampling rate of the signals ('dt').
    '''
    sigSum = None
    dtCheck = None
    for idx in range(len(varList)):
        sig, dt = extractStateVariable(mon, nIdx, varList[idx])
        if (idx == 0):
            sigSum = sig
            dtCheck = dt
        else:
            assert(dtCheck == dt)
            sigSum += sig

    return sigSum, dt



def extractSummedSignals(monList, varName, tStart, tEnd, monIdx=0):
    '''
    Extract state variables from a pair of monitors. One in the centre, the
    other one at the edge of the neural sheet.
    '''
    t, dt = extractStateVariable(monList, monIdx, 'times')
    sig, dt = sumAllVariables(monList, monIdx, varName)
    t, sig, idx = sliceSignal(t, sig, tStart, tEnd)
    return t, sig




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


class MonitoredTorusSpikes(spikes.TorusPopulationSpikes):
    '''Spikes that extract torus data from a NEST spike monitor.'''
    def __init__(self, data, monName, NXName, NYName):
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
        Nx = data['net_attr'][NXName]
        Ny = data['net_attr'][NYName]
        senders, times = extractSpikes(data[monName])
        super(MonitoredTorusSpikes, self).__init__(senders, times, (Nx, Ny))

