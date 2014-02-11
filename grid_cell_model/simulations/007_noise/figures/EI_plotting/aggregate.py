#
#   aggregate.py
#
#   Data aggregation, mainly parameter sweeps
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
from abc import ABCMeta, abstractmethod
import numpy as np
import numpy.ma as ma

from parameters import DataSpace
from otherpkg.log import log_warn

import logging
logger = logging.getLogger(__name__)


def aggregate2DTrial(sp, varList, trialNumList, fReduce=np.mean,
        ignoreNaNs=False):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space.

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    trialNumList : list of ints
        A list specifying exactly which trials are to be processed.
    fReduce : a function f(data, axis)
        A reduction function.
    ignoreNaNs : bool, optional
        If True, mask the NaN values.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=False)
    if (ignoreNaNs):
        nans = np.isnan(retVar)
        retVar = ma.MaskedArray(retVar, mask=nans)
    return fReduce(retVar, 2)


def aggregate2D(sp, varList, funReduce=None):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space, however the data is retrieved
    from the top-level of the data hierarchy, i.e. sp['analysis']. funReduce is
    applied on the data of each item in the ParamSpace (not necessarily
    trials).

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    funReduce : a function f(data, axis)
        A reduction function.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    return sp.aggregateData(varList, funReduce=funReduce,
            trialNumList='all-at-once', saveData=True)



def computeYX(sp, iterList, r=0, c=0, trialNum=0, normalize=True, **kw):
    E, I = sp.getIteratedParameters(iterList)
    if (normalize):
        Ne = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ne')
        Ni = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ni')
    else:
        Ne = 1
        Ni = 1
    return E/Ne, I/Ni



def computeVelYX(sp, iterList, r=0, c=0, trialNum=0, normalize=True, **kw):
    E, I = sp.getIteratedParameters(iterList)
    if (normalize):
        Ne = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0],
                'net_Ne')
        Ni = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0],
                'net_Ni')
    else:
        Ne = 1.0
        Ni = 1.0

    return E/Ne, I/Ni



def aggregateType(sp, iterList, types, NTrials, ignoreNaNs=False, **kw):
    '''
    Automatically aggregate data according to the type of the data.
    '''
    type, subType = types
    vars          = ['analysis']
    output_dtype  = 'array'
    funReduce     = None
    normalizeTicks = kw.pop('normalizeTicks', False)

    if (type == 'gamma'):
        # Gamma oscillation analyses
        if (subType == 'acVal'):
            # Autocorrelation first local maximum
            vars += ['acVal']
        elif (subType == 'freq'):
            # Gamm frequency
            vars += ['freq']
        elif (subType == 'acVec'):
            # All the autocorrelations
            vars += ['acVec']
            output_dtype = 'list'
        else:
            raise ValueError('Unknown gamma subtype: {0}'.format(subType))
        trialNumList  = np.arange(NTrials)

    elif (type == 'bump'):
        trialNumList  = np.arange(NTrials)
        if (subType == 'sigma'):
            vars += ['bump_e', 'sigma']
        elif (subType == 'rateMap_e'):
            vars += ['bump_e', 'bump_e_rateMap']
            output_dtype = 'list'
        elif (subType == 'rateMap_i'):
            vars += ['bump_i', 'bump_i_rateMap']
            output_dtype = 'list'
        else:
            raise ValueError('Unknown bump subtype: {0}'.format(subType))

    elif (type == 'bump_full'):
        trialNumList  = np.arange(NTrials)
        if (subType == 'sigma'):
            vars += ['bump_e_full', 'sigma']
        elif (subType == 'rateMap_e'):
            vars += ['bump_e_full', 'bump_e_rateMap']
            output_dtype = 'list'
        elif (subType == 'rateMap_i'):
            vars += ['bump_i_full', 'bump_i_rateMap']
            output_dtype = 'list'
        else:
            raise ValueError('Unknown bump_full subtype: {0}'.format(subType))

    elif (type == 'velocity'):
        if (subType == 'slope'):
            vars += ['lineFitSlope']
        elif (subType == 'fitErr'):
            vars += ['lineFitErr']
            funReduce = np.sum
        else:
            raise ValueError('Unknown velocity subtype: {0}'.format(subType))
        trialNumList = 'all-at-once'

    elif (type == 'grids'):
        trialNumList  = np.arange(NTrials)
        if (subType == 'gridnessScore'):
            vars += ['gridnessScore']
            funReduce = None
        else:
            raise ValueError('Unknown grids subtype: {0}'.format(subType))
    elif type == 'FR':
        trialNumList = np.arange(NTrials)
        if subType == 'E':
            vars += ['FR_e', 'avg']
            funReduce = None
        elif subType == 'I_10': # user should be aware of 10 neurons limitation
            vars += ['FR_i', 'all']
            funReduce = None
        else:
            raise ValueError('Unknown FR subtype: {0}'.format(subType))

    else:
        raise ValueError('Unknown aggregation type: {0}'.format(type))


    data = sp.aggregateData(vars, trialNumList, output_dtype=output_dtype,
            loadData=True, saveData=False, funReduce=funReduce)
    if (output_dtype != 'list'):
        data = ma.MaskedArray(data)
        if (ignoreNaNs):
            log_warn('aggregateType', 'Ignoring NaNs')
            nans = np.isnan(data)
            data.mask = nans

    if (type == 'velocity'):
        Y, X = computeVelYX(sp, iterList, normalize=normalizeTicks, **kw)
    else:
        if (type == 'bump' or type == 'bump_full'):
            if (subType == 'sigma'):
                # bump sigma is a reciprocal
                data = 1./data
                ignoreThreshold = 1.0
                data.mask = np.logical_or(data.mask, data > ignoreThreshold)
                data = np.mean(data, axis=2) # TODO: fix the trials, stack them
        else:
            data = np.mean(data, axis=2) # TODO: fix the trials, stack them
        Y, X = computeYX(sp, iterList, normalize=normalizeTicks, **kw)

    return data, X, Y



class AggregateData(object):
    __meta__ = ABCMeta
    analysisRoot = ['analysis']

    def __init__(self, space, iterList, NTrials, ignoreNaNs=False,
            normalizeTicks=False):
        self.sp = space
        self.iterList = iterList
        self.NTrials = NTrials
        self.ignoreNaNs = ignoreNaNs
        self.normalizeTicks = normalizeTicks


    @abstractmethod
    def getData(self):
        raise NotImplementedError()


bumpPosLogger = logging.getLogger('{0}.{1}'.format(__name__, 'BumpPositionData'))
class BumpPositionData(AggregateData):
    funReduce     = None
    output_dtype  = 'list'

    def __init__(self, space, iterList, NTrials, what, root, ignoreNaNs=False,
            normalizeTicks=False, **kw):
        super(BumpPositionData, self).__init__(space, iterList, NTrials,
                ignoreNaNs, normalizeTicks)
        self._what = what
        if self._what is None:
            raise ValueError('%s.what must be a list of strings' % self.what)
        self._root = root
        self._data = None
        self._vars = self.analysisRoot + self._root + self._what
        self._kw = kw
        bumpPosLogger.debug('self._vars: %s', self._vars)

    def getWhat(self): return self._what
    what = property(getWhat)

    def getData(self):
        if self._data is not None:
            return self._data

        trialNumList  = np.arange(self.NTrials)
        data = self.sp.aggregateData(self._vars,
                trialNumList,
                output_dtype=self.output_dtype,
                loadData=True,
                saveData=False,
                funReduce=None)
        Y, X = computeYX(self.sp, self.iterList, normalize=self.normalizeTicks,
                **self._kw)
        return data, X, Y


class AggregateBumpReciprocal(BumpPositionData):
    times_dtype = 'list'

    def __init__(self, space, iterList, NTrials, ignoreNaNs=False,
            normalizeTicks=True, root=['bump_e'], tStart=0, aggrFunc=np.median,
            **kw):
        what = ['positions', 'sigma']
        super(AggregateBumpReciprocal, self).__init__(space, iterList, NTrials,
                what, root, ignoreNaNs, normalizeTicks, **kw)
        self.tStart = tStart
        self.aggrFunc = aggrFunc

        trialNumList  = np.arange(self.NTrials)
        timeVars = self.analysisRoot + self._root + ['positions', 't']
        self._timeData = self.sp.aggregateData(timeVars,
                trialNumList,
                output_dtype=self.times_dtype,
                loadData=True,
                saveData=False,
                funReduce=None)[0][0][0]

    def getData(self):
        timeIdx = self._timeData >= self.tStart
        resShape = (self.sp.shape[0], self.sp.shape[1], self.NTrials)
        res = np.ma.MaskedArray(np.ndarray(resShape), mask=True)
        nRows, nCols = self.sp.shape

        rawData, X, Y = super(AggregateBumpReciprocal, self).getData()
        for r in xrange(nRows):
            for c in xrange(nCols):
                for trialIdx in xrange(self.NTrials):
                    reciprocal = 1./np.abs(rawData[r][c][trialIdx])
                    if isinstance(reciprocal, np.ndarray):
                        res[r, c, trialIdx] = \
                                self.aggrFunc(reciprocal[timeIdx])
                    elif not (isinstance(reciprocal, float) and
                            np.isnan(reciprocal)):
                        raise ValueError('Something went wrong here. ' + \
                                'reciprocal must be either array or NaN')
        return np.mean(res, axis=2), X, Y



def collapseSweeps(data):
    '''
    Take a list of 2D parameter sweep results, flatten all of them, and stack
    them vertically.
    '''
    stackedData = []
    for idx in xrange(len(data)):
        stackedData.append(data[idx].ravel())
    return np.ma.vstack(stackedData)


def collapseNoise(dataSpaces, iterList, types, NTrials, **kw):
    data = []
    for ns_idx, _ in enumerate(dataSpaces):
        d, X, Y = aggregateType(dataSpaces[ns_idx], iterList, types, NTrials,
                **kw)
        data.append(d)

    return collapseSweeps(data), X, Y
