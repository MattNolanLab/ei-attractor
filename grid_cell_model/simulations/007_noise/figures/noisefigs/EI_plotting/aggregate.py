'''
.. currentmodule:: EI_plotting.aggregate

Data aggregation, mainly parameter sweeps
'''
from abc import ABCMeta, abstractmethod
import numpy as np
import numpy.ma as ma

from grid_cell_model.parameters import DataSpace
from grid_cell_model.otherpkg.log import log_warn, getClassLogger
from grid_cell_model.analysis.image import Position2D
import grid_cell_model.analysis.image as image
import grid_cell_model.analysis.signal as asignal

import logging
logger = logging.getLogger(__name__)
gammaAggrLogger = getClassLogger('GammaAggregateData', __name__)


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
            normalizeTicks=False, collapseTrials=True, r=0, c=0):
        self.sp = space
        self.iterList = iterList
        self.NTrials = NTrials
        self.ignoreNaNs = ignoreNaNs
        self.normalizeTicks = normalizeTicks
        self.collapseTrials = collapseTrials
        self.r = r
        self.c = c

    @abstractmethod
    def getData(self):
        raise NotImplementedError()

    
    def filter_data(self, filter_obj):
        '''Apply a mask of a filter object to the current data.

        This operation returns a new AggregateData object
        '''
        return FilteredData(self, filter_obj)


class FilteredData(AggregateData):
    '''This applies a filter mask to a data object'''
    def __init__(self, data, filter_obj):
        self._data = data
        self._filter = filter_obj

    def getData(self):
        data, X, Y = self._data.getData()
        data = np.ma.array(data, copy=True)
        data.mask = np.logical_or(data.mask, self._filter.get_mask())
        return data, X, Y


class AggregateDataFilter(AggregateData):
    def __init__(self, data):
        self._data = data

    def getData(self):
        raise NotImplementedError()

    def filter_data(self, filter_obj):
        raise RuntimeError("Check if filter-on-filter operation makes sense")

    def get_mask(self):
        data, _, _ = self.getData()
        return np.copy(data.mask)


class NoZeroExcitationFilter(AggregateDataFilter):
    def __init__(self, data):
        super(NoZeroExcitationFilter, self).__init__(data)

    def getData(self):
        data, X, Y = self._data.getData()
        data = np.ma.array(data, copy=True)
        data.mask[:, 0] = True
        return data, X, Y


class GTFilter(AggregateDataFilter):
    '''Only use data that are greater than threshold'''
    def __init__(self, data, threshold):
        self.threshold = threshold
        super(GTFilter, self).__init__(data)

    def getData(self):
        data, X, Y = self._data.getData()
        data = np.ma.array(data, copy=True)
        data.mask = np.logical_or(data.mask, data <= self.threshold)
        return data, X, Y


def maskNaNs(a, really):
    return np.ma.masked_invalid(a, copy=False) if really else a


##############################################################################
# Grids / score
class GridnessScore(AggregateData):
    '''Extract gridness score from the aggregated data'''
    def __init__(self, space, iterList, **kw):
       super(GridnessScore, self).__init__(space, iterList, None, **kw)
       self._gscore = None

    def _getRawData(self):
        if self._gscore is None:
            path = self.analysisRoot + ['gridnessScore']
            self._gscore = self.sp.getReduction(path)
            self._Y, self._X = computeYX(self.sp, self.iterList,
                                         normalize=self.normalizeTicks,
                                         r=self.r, c=self.c)
        return self._gscore, self._X, self._Y

    def getData(self):
        data, X, Y = self._getRawData()
        return np.mean(maskNaNs(data, self.ignoreNaNs), axis=2), X, Y

        
class GammaAggregateData(AggregateData):
    '''Extract power of gamma oscillations from the aggregated data'''
    def __init__(self, what, space, iterList, **kw):
        '''``what`` determines the data field to extract'''
        super(GammaAggregateData, self).__init__(space, iterList, None, **kw)
        self._acval = None
        self._what = what

    def _getRawData(self):
        if self._acval is None:
            path = self.analysisRoot + [self._what]
            gammaAggrLogger.info("Extracting data from path: %s", path)
            self._acval = self.sp.getReduction(path)
            self._Y, self._X = computeYX(self.sp, self.iterList,
                                         normalize=self.normalizeTicks)
        return self._acval, self._X, self._Y

    def getData(self):
        data, X, Y = self._getRawData()
        if self.collapseTrials:
            return np.mean(maskNaNs(data, self.ignoreNaNs), axis=2), X, Y
        else:
            return data, X, Y



##############################################################################
# Bumps
class IsBump(AggregateData):
    '''Retrieve bump classification data from the space.'''
    
    def __init__(self, space, iterList, ignoreNaNs=False, normalizeTicks=False,
                 **kw):
        super(IsBump, self).__init__(space, iterList, None, ignoreNaNs,
                normalizeTicks, **kw)
        self._fracTotal = None
        self._X = None
        self._Y = None

    def _getRawData(self):
        if self._fracTotal is None:
            path = self.analysisRoot + ['bump_e/isBump/fracTotal']
            self._fracTotal = self.sp.getReduction(path)
            self._Y, self._X = computeYX(self.sp, self.iterList,
                    self.normalizeTicks)
        return self._fracTotal, self._X, self._Y

    def getData(self):
        data, X, Y = self._getRawData()
        if self.collapseTrials:
            return np.mean(maskNaNs(data, self.ignoreNaNs), axis=2), X, Y
        else:
            return data, X, Y


class IsBumpCollapsed(AggregateData):
    '''Take a list of data spaces and aggregate them and their trials (fracTotal
    data set) into one vector usable for histogram calculation.'''

    def __init__(self, spaces, iterList, ignoreNaNs=False,
            normalizeTicks=False):
        self.spaces = spaces
        super(IsBumpCollapsed, self).__init__(None, iterList, None,
                ignoreNaNs, normalizeTicks)
        self.fracList = []
        for sp in spaces:
            self.fracList.append(IsBump(sp, iterList, ignoreNaNs,
                normalizeTicks))


    def getData(self):
        collapsedData = []
        for fracObj in self.fracList:
            fracTotal, X, Y = fracObj.getData()
            collapsedData.append(fracTotal.ravel())
        return np.hstack(collapsedData), X, Y



class BumpFormationFilter(IsBump):
    '''Provide a filter that  will mask data where bump is not formed.

    The threshold of the masking process must be supplied.
    '''
    def __init__(self, threshold, *args, **kwargs):
        super(BumpFormationFilter, self).__init__(*args, **kwargs)
        self.threshold = threshold
        self._filter = None
        self._X     = None
        self._Y     = None

    def getData(self):
        if self._filter is None:
            fracTotal, self._X, self._Y = \
                    super(BumpFormationFilter, self).getData()
        self._filter = fracTotal > self.threshold
        return self._filter, self._X, self._Y

    def filterData(self, data, invert=False):
        '''
        if `invert`, then return data which do *not* match the filter.
        '''
        filter, _, _ = self.getData()
        if filter.shape != data.shape:
            raise TypeError("Shapes of data to filter and filter must match.")
        filt = np.logical_not(filter) if invert else filter
        mask = np.logical_or(getattr(data, 'mask', False), filt)
        return np.ma.MaskedArray(data, mask=mask)
        





bumpPosLogger = logging.getLogger('{0}.{1}'.format(__name__, 'BumpPositionData'))
class BumpPositionData(AggregateData):
    funReduce     = None
    output_dtype  = 'list'
    times_dtype = 'list'

    def __init__(self, space, iterList, NTrials, what, root, ignoreNaNs=False,
            normalizeTicks=False, **kw):
        super(BumpPositionData, self).__init__(space, iterList, NTrials,
                ignoreNaNs, normalizeTicks)
        self._what = what
        if self._what is None:
            self._what = []
        self._root = root
        self._data = None
        self._X = None
        self._Y = None
        self._vars = self.analysisRoot + self._root + self._what
        self._kw = kw
        bumpPosLogger.debug('self._vars: %s', self._vars)

        trialNumList  = np.arange(self.NTrials)
        timeVars = self.analysisRoot + self._root + ['positions', 'times']
        self._timeData = self.sp.aggregateData(timeVars,
                trialNumList,
                output_dtype=self.times_dtype,
                loadData=True,
                saveData=False,
                funReduce=None)[0][0][0]

    def getWhat(self): return self._what
    what = property(getWhat)


    def _getRawData(self, vars):
        trialNumList  = np.arange(self.NTrials)
        data = self.sp.aggregateData(vars,
                trialNumList,
                output_dtype=self.output_dtype,
                loadData=True,
                saveData=False,
                funReduce=None)
        if self._X is None: # or self._Y is None
            self._Y, self._X = computeYX(self.sp, self.iterList,
                    normalize=self.normalizeTicks, **self._kw)
        return data, self._X, self._Y


    def getData(self):
        if self._data is None:
            self._data, self._X, self._Y = self._getRawData(self._vars)
        return self._data, self._X, self._Y

    def getTimes(self):
        return self._timeData


bumpDiffLogger = logging.getLogger('{0}.{1}'.format(__name__,
        'BumpDifferencePosition'))
class BumpDifferencePosition(BumpPositionData):
    '''
    Compute vectors of distances from the snapshot of the bump specified by a
    user argument. Everything before the start time will be ignored.
    '''
    def __init__(self, space, iterList, NTrials, ignoreNaNs=False,
            normalizeTicks=True, root=['bump_e'], tStart=0, **kw):
        what = None # This should not be valid
        super(BumpDifferencePosition, self).__init__(space, iterList, NTrials,
                what, root, ignoreNaNs, normalizeTicks, **kw)
        self.tStart = tStart
        bumpDiffLogger.warn('Nx, Ny are fixed in the code. Make sure the '+\
                "torus size is the same as specified here.")


    def _getMus(self):
        varBase = self.analysisRoot + self._root + ['positions']
        mu_x_all, X, Y = self._getRawData(varBase + ['mu_x'])
        mu_y_all, _, _ = self._getRawData(varBase + ['mu_y'])
        return mu_x_all, mu_y_all, X, Y


    def getData(self):
        mu_x_all, mu_y_all, X, Y = self._getMus()
        timeIdx = self._timeData >= self.tStart
        nRows, nCols = self.sp.shape
        distances = [[[None for trialIdx in xrange(self.NTrials)] for c in
                xrange(nCols)] for r in xrange(nRows)]

        for r in xrange(nRows):
            for c in xrange(nCols):
                #if r == 13 and c == 20:
                #    import pdb; pdb.set_trace()
                for trialIdx in xrange(self.NTrials):
                    try:
                        mu_x = mu_x_all[r][c][trialIdx]
                        mu_y = mu_y_all[r][c][trialIdx]
                    except KeyError:
                        import pdb; pdb.set_trace()
                    if isinstance(mu_x, np.ndarray) and isinstance(mu_y,
                            np.ndarray):
                        mu_x = mu_x[timeIdx]
                        mu_y = mu_y[timeIdx]

                        startPos = Position2D(mu_x[0], mu_y[0])
                        positions = Position2D(mu_x, mu_y)
                        torusSize = Position2D(34, 30)
                        distances[r][c][trialIdx] = image.remapTwistedTorus(\
                                startPos, positions, torusSize)

        return distances, X, Y
                    

    def getTimes(self):
        return self._timeData[self._timeData >= self.tStart]


class BumpAvgDifferenceFromPos(BumpDifferencePosition):
    def __init__(self, startPos, *args, **kwargs):
        super(BumpAvgDifferenceFromPos, self).__init__(*args, **kwargs)
        self.startPos = startPos
        self._avgDiff = None
        self._X, self._Y = None, None


    def getData(self):
        if self._avgDiff is not None:
            return self._avgDiff, self._X, self._Y

        distances, self._X, self._Y = super(BumpAvgDifferenceFromPos,
                self).getData()
        self._avgDiff = np.ma.MaskedArray(np.ndarray(self.sp.shape), mask=True)
        nRows, nCols = self.sp.shape
        for r in xrange(nRows):
            for c in xrange(nCols):
                trialDiffs = []
                for trialIdx in xrange(self.NTrials):
                    d = distances[r][c][trialIdx]
                    if d is not None and isinstance(d, np.ndarray):
                        mn = np.mean(d)
                        if not np.isnan(mn):
                            trialDiffs.append(mn)
                if len(trialDiffs) > 0:
                    self._avgDiff[r, c] = np.mean(trialDiffs)
        return self._avgDiff, self._X, self._Y



class BumpDifferenceAtTime(BumpDifferencePosition):
    def __init__(self, startPos, diffTime, *args, **kwargs):
        super(BumpDifferenceAtTime, self).__init__(*args, tStart=0, **kwargs)
        self.startPos = startPos
        self.diffTime = diffTime
        bumpDiffLogger.warn('Nx, Ny are fixed in the code. Make ' + \
                "sure the torus size is the same as specified here.")

    def getData(self):
        diffIdx = np.nonzero(self.getTimes() >= self.diffTime)[0]
        if len(diffIdx) == 0:
            msg = 'Cannot find appropriate index in the bump position data. '+\
                    'Check your difference time: {0}'
            raise IndexError(msg.format(self.diffTime))
        #diffIdx = diffIdx[-1]

        mu_x_all, mu_y_all, X, Y = self._getMus()
        nRows, nCols = self.sp.shape
        diffs = np.ma.MaskedArray(np.ndarray(self.sp.shape), mask=True)
        for r in xrange(nRows):
            for c in xrange(nCols):
                trialDiffs = []
                for trialIdx in xrange(self.NTrials):
                    mu_x = mu_x_all[r][c][trialIdx]
                    mu_y = mu_y_all[r][c][trialIdx]
                    if isinstance(mu_x, np.ndarray) and isinstance(mu_y,
                            np.ndarray):
                        pos_x = mu_x[diffIdx]
                        pos_y = mu_y[diffIdx]
                        diffAtTLogger.warn('Nx, Ny are fixed in the code. Make ' + \
                                "sure the torus size is the same as specified here.")
                        torusSize = Position2D(34, 30)
                        d = image.remapTwistedTorus(\
                                Position2D(self.startPos[0], self.startPos[1]),
                                Position2D(np.array([pos_x]), np.array([pos_y])),
                                torusSize)
                        trialDiffs.append(d[0])
                if len(trialDiffs) > 0:
                    diffs[r, c] = np.mean(trialDiffs) 
        return diffs, X, Y

diffAtTLogger = logging.getLogger("{0}.{1}".format(__name__,
        'BumpDifferenceAtTime'))




class BumpDriftAtTime(BumpDifferencePosition):
    def __init__(self, tDrift, *args, **kwargs):
        super(BumpDriftAtTime, self).__init__(*args, **kwargs)
        self.tDrift = tDrift


    def getData(self):
        # Find the correct bin if it exists. Time array might not start with
        # t=0!
        driftIdx = np.nonzero(self.getTimes() <= self.tDrift)[0]
        if len(driftIdx) == 0:
            msg = 'Cannot find appropriate index in the bump position data. '+\
                    'Check your drift time: {0}'
            raise IndexError(msg.format(self.tDrift))
        driftIdx = driftIdx[-1]

        drifts = np.ma.MaskedArray(np.ndarray(self.sp.shape), mask=True)
        distances, X, Y = super(BumpDriftAtTime, self).getData()
        nRows, nCols = self.sp.shape
        for r in xrange(nRows):
            for c in xrange(nCols):
                trialDst = []
                for trialIdx in xrange(self.NTrials):
                    dist = distances[r][c][trialIdx]
                    if dist is not None:
                        trialDst.append(dist[driftIdx])
                if len(trialDst) != 0:
                    drifts[r, c] = np.mean(trialDst)
                    if drifts[r, c] > 20:
                        msg = "drift > 20: {0}, r: {1}, c:{2}"
                        diffAtTLogger.info(msg.format(drifts[r, c], r, c))
        #totalTime = (self.tDrift - self.tStart) * 1e-3
        return drifts, X, Y
                    



class AggregateBumpReciprocal(BumpPositionData):

    def __init__(self, space, iterList, NTrials, ignoreNaNs=False,
            normalizeTicks=True, root=['bump_e'], tStart=0, aggrFunc=np.median,
            **kw):
        what = ['positions', 'sigma']
        super(AggregateBumpReciprocal, self).__init__(space, iterList, NTrials,
                what, root, ignoreNaNs, normalizeTicks, **kw)
        self.tStart = tStart
        self.aggrFunc = aggrFunc
        self._trialData = None
        self._X = None
        self._Y = None

    def _computeData(self):
        if self._trialData is None:
            timeIdx = self._timeData >= self.tStart
            resShape = (self.sp.shape[0], self.sp.shape[1], self.NTrials)
            self._trialData = np.ma.MaskedArray(np.ndarray(resShape), mask=True)
            nRows, nCols = self.sp.shape

            rawData, self._X, self._Y = super(AggregateBumpReciprocal, self).getData()
            for r in xrange(nRows):
                for c in xrange(nCols):
                    for trialIdx in xrange(self.NTrials):
                        reciprocal = 1./np.abs(rawData[r][c][trialIdx])
                        if isinstance(reciprocal, np.ndarray):
                            self._trialData[r, c, trialIdx] = \
                                    self.aggrFunc(reciprocal[timeIdx])
                        elif not (isinstance(reciprocal, float) and
                                np.isnan(reciprocal)):
                            raise ValueError('Something went wrong here. ' + \
                                    'reciprocal must be either array or NaN')
            self._trialData

        return self._trialData, self._X, self._Y


    def getData(self):
        trialData, X, Y = self._computeData()
        return np.mean(trialData, axis=2), X, Y

    def getTrialData(self):
        trialData, _, _ = self._computeData()
        return trialData

    def getTimes(self):
        return self._timeData[self._timeData >= self.tStart]




##############################################################################
#                            Firing rates


class PopulationFR(AggregateData):

    def __init__(self, space, iterList, **kw):
        super(PopulationFR, self).__init__(space, iterList, None, **kw)
        self._FR = None
        self._X  = None
        self._Y  = None

    def _getRawData(self):
        if self._FR is None:
            path = self.analysisRoot[0] + '/FR_e/popSliding'
            FR = self.sp.getReduction(path)
            nTrials = len(FR[0][0])
            dataLen = len(FR[0][0][0])
            self._FR = np.ndarray((self.sp.shape[0], self.sp.shape[1],
                                        nTrials, dataLen),
                                     dtype=np.double)
            for r in xrange(self.sp.shape[0]):
                for c in xrange(self.sp.shape[1]):
                    for trialNum in xrange(nTrials):
                        data = FR[r][c][trialNum]
                        self._FR[r, c, trialNum, :] = \
                                data if data is not None else np.nan

            self._Y, self._X = computeYX(self.sp, self.iterList,
                                         normalizeTicks=self.normalizeTicks)

        return self._FR, self._X, self._Y


class MaxPopulationFR(PopulationFR):
    '''
    Extract the maximal (sliding) population firing rate for each parameter
    setting.
    '''

    def __init__(self, space, iterList, **kw):
        super(MaxPopulationFR, self).__init__(space, iterList, **kw)

    def getData(self):
        FR, X, Y = self._getRawData()
        maxFR = np.nanmax(FR, axis=3)
        maxFR[maxFR == 5000] = np.nan
        return (np.mean(maskNaNs(maxFR, self.ignoreNaNs), axis=2), self._X,
                self._Y)



class MaxThetaPopulationFR(PopulationFR):
    '''
    Extract the median/mean of the maximal population firing rate *every theta*
    cycle.
    '''
    def __init__(self, thetaT, sig_dt, reduceFunc, space, iterList, **kw):
        super(MaxThetaPopulationFR, self).__init__(space, iterList, **kw)
        self.thetaT = thetaT
        self.sig_dt = sig_dt
        self._trialMax = None        # Trial data, non-reduced
        self.reduceFunc = reduceFunc

    def getNonReducedData(self):
        if self._trialMax is None:
            FR, _, _ = self._getRawData()
            nTrials = FR.shape[2]
            self._trialMax = np.ndarray((self.sp.shape[0], self.sp.shape[1],
                                            nTrials),
                                         dtype=object)
            for r in xrange(self.sp.shape[0]):
                for c in xrange(self.sp.shape[1]):
                    for trialNum in xrange(nTrials):
                        rate = FR[r, c, trialNum, :]
                        self._trialMax[r, c, trialNum] = np.max(
                                asignal.splitSigToThetaCycles(rate,
                                                              self.thetaT,
                                                              self.sig_dt),
                                axis=1)
        return self._trialMax

    def getData(self):
        data = self.getNonReducedData()
        reduction = np.vectorize(self.reduceFunc, cache=True)
        reducedData = reduction(data)
        return (np.mean(maskNaNs(reducedData, self.ignoreNaNs), axis=2),
                self._X, self._Y)


class AvgPopulationFR(AggregateData):
    '''Average population firing rate, averaged over all neurons in the
    network'''
    def __init__(self, where, space, iterList, **kw):
        super(AvgPopulationFR, self).__init__(space, iterList, None, **kw)
        self._where = where
        self._FR = None
        self._X  = None
        self._Y  = None

    def _getRawData(self):
        if self._FR is None:
            path = "{root}/{where}".format(
                        root=self.analysisRoot[0],
                        where = self._where)
            self._FR = np.asarray(self.sp.getReduction(path))
            self._Y, self._X = computeYX(self.sp, self.iterList,
                                         normalizeTicks=self.normalizeTicks,
                                         r=self.r, c=self.c)

        return self._FR, self._X, self._Y

    def getData(self):
        data, X, Y = self._getRawData()
        return np.mean(maskNaNs(data, self.ignoreNaNs), axis=2), X, Y


##############################################################################
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


def collapseNoiseAggregated(dataList):
    data = []
    for dataObject in dataList:
        d, X, Y = dataObject.getData()
        data.append(d)

    return collapseSweeps(data), X, Y

