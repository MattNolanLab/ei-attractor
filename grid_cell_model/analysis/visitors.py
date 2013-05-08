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
import numpy as np
from otherpkg.log     import log_info
from analysis.signal  import localExtrema, butterBandPass, autoCorrelation

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

    def visitDictDataSet(self, ds):
        '''
        Visit the dictionary data set, 'ds', and perform the specific operations
        (defined by the derived classes) on this data set.
        '''
        raise NotImplementedError()




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



def findFreq(ac, dt, ext_idx, ext_t):
    '''
    Find the first local maximum in an autocorrelation function and extract the
    frequency and the value of the autocorrelation at the detected frequency.

    Parameters
    ----------
    ac : numpy vector
        A vector containing the autocorrelation function
    dt : float
        Sampling rate
    ext_idx : numpy array
        An array containig the indexes of local extrema (both maxima and
        minima) in 'ac'
    ext_t : numpy array
        The array of the same size as 'ext_idx', which contains the types of
        the local extrema ( >0 ... max, <0 ... min)
    output
        A tuple containig ('freq', 'corr'), where 'freq' is the frequency at
        the first local maximum and 'corr' is the value of the autocorrelation
        function at the same point.
    '''
    max_idx = np.nonzero(ext_t > 0)[0]
    if (len(max_idx) == 0):
        return (np.nan, np.nan)
    
    # First local maximum ac[0] excluded
    max1_idx = ext_idx[max_idx[0]]
    max1_t   = max1_idx * dt
    max1     = ac[max1_idx]

    return (1./max1_t, max1)




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



class AutoCorrelationVisitor(DictDSVisitor):
    '''
    A visitor to compute autocorrelations of state monitor data and extract
    information from them.

    The autocorrelation visitor takes as an input a state monitor, computes
    autocorrelations (with specified lag) of the specified synaptic currents
    and detects the major frequency, and power at that frequency, for all the
    monitored neurons. The results will be stored to the dictionary where the
    data came from.
    '''
    def __init__(self, monName, stateList, dtMult=1e-3, norm=True,
            bandStart=20, bandEnd=200, forceUpdate=False):
        '''
        Initialise the visitor.

        Parameters
        ----------
        monName : string
            Name of the monitor; key in the data set dictionary
        stateList : list of strings
            A list of strings naming the state variables to extract (and sum)
        dtMult : float, optional
            dt Multiplier to transform dt into seconds
        norm : bool, optional
            Whether the autocorrelation function should be normalized
        bandStart : float, optional
            Bandpass start frequency
        bandEnd   : float, optional
            Bandpass end frequency
        forceUpdate : bool
            Whether to compute and store all the data even if they already
            exist in the data set.
        '''
        self.monName     = monName
        self.stateList   = stateList
        self.maxLag      = None
        self.dtMult      = dtMult
        self.norm        = norm
        self.bandStart   = bandStart
        self.bandEnd     = bandEnd
        self.forceUpdate = forceUpdate



    def extractACStat(self, mon):
        '''
        Extract autocorrelation statistics from a monitor.

        For each monitored neuron, extract the (highest) frequency, value of
        the autocorrelation at the frequency and the autocorrelation function
        itself.
    
        Parameters
        ----------
        mon : list of dicts
            A list of (NEST) state monitors' status dictionaries
        output : tuple
            A tuple (freq, acval, acVec), containing the arrays of frequencies
            for the monitored neurons, autocorrelation values at the
            corresponding frequencies, and autocorrelation functions of all the
            neurons.
        '''
        freq   = [] # Frequency of input signal
        acval  = [] # Auto-correlation at the corresponding frequency
        acVec  = []
        for n_id in range(len(mon)):
        #for n_id in range(5):
            #print "n_id: ", n_id
            sig, dt = sumAllVariables(mon, n_id, self.stateList)
            sig = butterBandPass(sig, dt*self.dtMult, self.bandStart,
                    self.bandEnd)
            ac = autoCorrelation(sig - np.mean(sig), max_lag=self.maxLag/dt,
                    norm=self.norm)
            ext_idx, ext_t = localExtrema(ac)
            acVec.append(ac)
    
            f, a = findFreq(ac, dt*1e-3, ext_idx, ext_t)
            freq.append(f)
            acval.append(a)
    
        return freq, acval, acVec


    def visitDictDataSet(self, ds):
        '''
        Visit the dictionary data set and extract frequency, autocorrelation
        for the detected frequency, and autocorrelation functions, for all the
        monitored neurons.  The parameters are defined by the constructor of
        the object.

        If the analysed data is already present, the analysis and storage of
        the data will be skipped.

        Parameters
        ----------
        ds : a dict-like object
            A data set to perform analysis on.
        '''
        data = ds.data
        if (('analysis' not in data.keys()) or self.forceUpdate):
            log_info("visitors", "Analysing a dataset")
            o = data['options']
            self.maxLag = 1. / (o['theta_freq'] * 1e-3)
            freq, acVal, acVec = self.extractACStat(data[self.monName])
            data['analysis'] = {
                    'freq'  : np.array(freq),
                    'acVal' : np.array(acVal),
                    'acVec' : np.array(acVec)
            }
        else:
            log_info("visitors", "Data present. Skipping analysis.")

