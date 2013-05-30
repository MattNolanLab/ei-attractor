#
#   visitors.py
#
#   Data analysis visitors. 
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

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



def extractSpikes(mon):
    '''
    Extract spikes from a spike monitor (a dict-like object), that contains the
    relevant fields.
    
    Return a tuple (senders, spikeTimes).
    '''
    e = mon['events']
    return (e['senders'], e['times'])


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





class Visitor(object):
    '''
    An abstract visitor class.

    Normally, the base class of the Visitor design pattern must contain all the
    methods implemented in the derived classes. Due to duck typing, one does
    not need to declare the specific implementation methods of the visitor.
    '''
    def __init__(self):
        raise NotImplementedError()



class DictDSVisitor(Visitor):
    '''
    Dictionary data set visitor.

    A visitor that takes a dictionary data set method of any kind and processes
    it. All the keys in the dictionary must be strings.
    '''
    def __init__(self):
        raise NotImplementedError()

    def visitDictDataSet(self, ds, **kw):
        '''
        Visit the dictionary data set, 'ds', and perform the specific operations
        (defined by the derived classes) on this data set.
        '''
        raise NotImplementedError()

    def getOption(self, data, optStr):
        '''Extract an option from a data dictionary'''
        return data['options'][optStr]

    def getNetParam(self, data, p):
        '''Extract a network parameter (p) from the data dictionary'''
        return data['net_attr'][p]


    def folderExists(self, d, nameList):
        '''
        Check if the folder at the end of the 'nameList' exists in dictionary
        'd'
        '''
        for name in nameList:
            if (name in d.keys()):
                d = d[name]
            else:
                return False
        return True

    def _checkAttrIsNone(self, attr, pName, data):
        '''
        Check if 'attr' is None. If yes, extract it from 'data' under the name
        'pName', otherwise return attr back.
        '''
        if (attr is None):
            return self.getOption(data, pName)
        else:
            return attr


    def _getSpikeTrain(self, data, monName, dimList):
        '''
        Return the senders and spike times from a monitor in the data
        dictionary

        Parameters
        ----------
        data : dict
            A data dictionary containing monName
        monName : str
            The name of the monitor
        dimList : list of string
            A  list of dimensions. This will be used to extract the total
            number of neurons from network parameters.
        output: tuple
            A tuple containing (senders, times, N). 'N' is the total number of
            neurons in the population.
        '''
        N = 1
        for dim in dimList:
            dimN = self.getNetParam(data, dim)
            N *= dimN
        mon = data[monName]
        senders, times = extractSpikes(mon)
        return senders, times, N


