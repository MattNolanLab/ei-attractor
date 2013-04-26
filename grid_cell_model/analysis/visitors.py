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
    '''Abstract visitor class'''
    def __init__(self):
        raise NotImplementedError()



class DictDSVisitor(Visitor):
    '''Dictionary data set visitor'''
    def __init__(self):
        raise NotImplementedError()

    def visitDictDataSet(self, ds):
        raise NotImplementedError()




def extractStateVariable(mon, nIdx, varStr):
    '''Extract state variable from a monitor.
    
    Parameters
    ----------
    mon : dict
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
    max_idx = np.nonzero(ext_t > 0)[0]
    if (len(max_idx) == 0):
        return (np.nan, np.nan)
        #raise ValueError("Autocorrelation must contain at least one local maximum")
    
    # First local maximum ac[0] excluded
    max1_idx = ext_idx[max_idx[0]]
    max1_t   = max1_idx * dt
    max1     = ac[max1_idx]

    return (1./max1_t, max1)




def sumAllVariables(mon, nIdx, varList):
    '''
    Extract all variables from the list of monitors and sum them
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
    def __init__(self, monName, stateList, dtMult=1e-3, norm=True,
            bandStart=20, bandEnd=200, forceUpdate=False):
        self.monName     = monName
        self.stateList   = stateList
        self.maxLag      = None
        self.dtMult      = dtMult
        self.norm        = norm
        self.bandStart   = bandStart
        self.bandEnd     = bandEnd
        self.forceUpdate = forceUpdate


    def visitDictDataSet(self, ds):
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

    def extractACStat(self, mon):
        '''
        Extrac autocorrelation statistics from a monitor
    
        Parameters
        ----------
        mon : list of dicts
            A list of (NEST) state monitors' status dictionary
        stateList : list of strings
            A list of strings naming the state variables to extract (and sum)
        maxLag : float
            Maximal lag to extract (in seconds, NOT timesteps)
        dtMult : float, optional
            dt Multiplier to transform dt into seconds
        norm : bool, optional
            Whether the autocorrelation function should be normalized
        bandStart : float, optional
            Bandpass start frequency
        bandEnd   : float, optional
            Bandpass end frequency
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
    
            #plt.figure()
            #plt.plot(times[0:slice], ac[0:slice])
            #plt.hold('on')
            #plt.plot(times[ext_idx], ac[ext_idx], '.')
            #plt.ylim([-1, 1])
            #plt.savefig(output_fname + "_peaks_ac_extrema_%d.pdf" % n_id)
            #plt.close()
    
            f, a = findFreq(ac, dt*1e-3, ext_idx, ext_t)
            freq.append(f)
            acval.append(a)
    
        return freq, acval, acVec


