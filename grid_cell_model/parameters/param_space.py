#
#   param_space.py
#
#   Parameter spaces classes.
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
from collections    import Sequence
from otherpkg.log   import log_warn, log_info

from data_storage       import DataStorage
from data_storage.dict  import getDictData
from data_sets          import DictDataSet

__all__ = []



class DataSpace(Sequence):
    '''
    An interface to the space of data with different parameters.
    '''

    def __init__(self, vals, pList=None, key=None):
        self._key  = key
        self._vals = vals
        self._pList = pList


    def __getitem__(self, key):
        return self._vals[key]


    def visit(self, visitor):
        for val in self._vals:
            val.visit(visitor)

    def __len__(self):
        return len(self._vals)

    @staticmethod 
    def getParam(data, paramStr):
        return data['options'][paramStr]

    @staticmethod
    def getNetParam(data, p):
        return data['net_attr'][p]



class TrialSet(DataSpace):
    '''
    A 1D DataSpace that contains a list of DataSet objects.
    '''

    def __init__(self, fileName, fileMode):
        self._fileName = fileName
        try:
            self._ds = DataStorage.open(fileName, fileMode)
            DataSpace.__init__(self, self._ds['trials'], key='trials')
            log_info("param_space", "Opened " + fileName)
        except IOError as e:
            self._ds = None
            msg =  "Could not open file {0}. Creating an empty DataSet instead."
            log_warn("param_space", msg.format(fileName))
            log_warn("param_space", "IOError message: {0}".format(str(e)))
            DataSpace.__init__(self, [], key='trials')

    def __del__(self):
        if (self._ds is not None):
            log_info("param_space", "Closing: " + self._fileName)
            self._ds.close()

    def __getitem__(self, key):
        return DictDataSet(self._vals[key])

    def getAllTrialsAsDataSet(self):
        return DictDataSet(self._ds)

    def visit(self, visitor, trialList=None):
        '''
        Apply a visitor to the trials set. There are two modes: Apply the
        visitor to trials separately, or pass the raw DictDataset to the
        visitor so it can apply some reduction operation onto the data set. See
        the description of the trialList parameter.

        Parameters
        ----------
        visitor : a DictDataSetVisitor object
            The visitor to apply on this data space.
        trialList : list, or None, or string
            If trialList is None, apply the visitor to all the trials
            separately. If trialList == 'all-at-once', pass on all the trials
            in one DataSet object to the visitor. Otherwise trialList must be a
            list of trials to pass on to the visitor (separately).
        '''
        if (trialList == 'all-at-once'):
            self.getAllTrialsAsDataSet().visit(visitor, fileName=self._fileName)
        else:
            if (trialList is None):
                trialList = xrange(len(self))

            for trialIdx in trialList:
                trial = self[trialIdx]
                trial.visit(visitor, fileName=self._fileName, trialNum=trialIdx)


class DummyTrialSet(DataSpace):
    '''An empty DataSpace that does nothing.'''
    def __init__(self):
        DataSpace.__init__(self, [], key='trials')

    def __getitem__(self, key):
        raise ValueError("A DummyTrialSet cannot be indexed!")

    def visit(self, visitor, trialList=None):
        pass



class JobTrialSpace2D(DataSpace):
    def __init__(self, shape, rootDir, dataPoints=None, fileMode='r+',
            fileFormat="job{0:05}_output.h5", forceWMode=False,
            checkParams=False):
        if (fileMode == 'w' and forceWMode == False):
            raise ValueError("'w' file open mode is not allowed. Use " +
                    "'forceWMode' to override.")
        self._shape = shape
        self._rootDir = rootDir
        self._dataPoints = dataPoints
        if (self._dataPoints is not None):
            self._partial = True
        else:
            self._partial = False
        self._fileFormat = fileFormat
        self._checkParams = checkParams
        self.rows = shape[0]
        self.cols = shape[1]
        rowData = []
        it = 0
        for row in xrange(self.rows):
            colData = []
            for col in xrange(self.cols):
                # Here either we have a full data space, or a particular row
                # and column are present in the list of selected data points.
                # Otherwise, append an empty list, i.e. this will do nothing
                if (self._dataPoints is None or
                        (row, col) in self._dataPoints):
                    fileName = rootDir + '/' + fileFormat.format(it)
                    colData.append(TrialSet(fileName, fileMode))
                else:
                    colData.append(DummyTrialSet())
                it += 1
            rowData.append(DataSpace(colData))
        DataSpace.__init__(self, rowData)
            
                
    def __len__(self):
        return self._shape[0] 


    def getShape(self):
        return self._shape


    def getIteratedParameters(self, nameList):
        if (len(nameList) != 2):
            raise ValueError("nameList must contain exactly 2 elements.")
        iterFileName = "{0}/iterparams.h5".format(self._rootDir)
        ds = DataStorage.open(iterFileName, 'r')
        ret = []
        for nm in nameList:
            ret.append(np.reshape(ds['iterParams'][nm], self._shape))
            if (self._checkParams):
                self._checkIteratedParameters(nm, ret[-1])
        ds.close()
        return ret


    def _checkIteratedParameters(self, paramStr, toCheck):
        tol  = 1e-9 * np.min(toCheck.flatten())
        msgStr = "Parameter {0}:[{1}][{2}: {3}] does not match."
        for r in xrange(self.rows):
            for c in xrange(self.cols):
                trials = self[r][c]
                for trial_idx in xrange(len(trials)):
                    pVal = self.getParam(trials[trial_idx].data, paramStr)
                    err = np.abs(pVal - toCheck[r][c])
                    if (err > tol):
                        raise Exception(msgStr.format(paramStr, r, c, err))

    def visit(self, visitor, trialList=None):
        for r in xrange(self.rows):
            for c in xrange(self.cols):
                self[r][c].visit(visitor, trialList)


    def aggregateData(self, varList, trialNumList, funReduce=None):
        '''
        Aggregate the data from each trial into a 2D object array of the shape
        (row, col), each item containing a list of values, one value for each
        trial.

        Parameters
        ----------
        do : dict-like object, or None
            Dictionary output. If None, return the result only
        name : string
            A key under which to store the object. If 'do' is None, this can be
            None as well
        varList : list of strings
            A path to the hierarchical dictionary
        trialNumList : list of ints or string
            A list of trials to aggregate. If string and the value is
            'all-at-once', aggregate the variable from the top level hierarchy
            of data.
        funReduce : a function f(x), or None
            A function to apply to each data point for each trial. Must take
            exactly one parameter. If None, no function will be applied.
        output : A 3D numpy array if trialNumList is a list, or a 2D array
                 otherwise
            All the aggregated data
        '''
        if (self._partial):
            # Cannot do aggregation on a restricted data set
            raise NotImplementedError("Data aggregation on a partial data " +
                    "space has not been implemented yet.")
        shape = self.getShape()
        rows = shape[0]
        cols = shape[1]
        retVar = None
        if (trialNumList == 'all-at-once'):
            retVar = np.ndarray((rows, cols))
        else:
            nTrials = len(trialNumList)
            retVar = np.ndarray((rows, cols, nTrials))

        if (funReduce is None):
            funReduce = lambda x: x
        for r in xrange(rows):
            for c in xrange(cols):
                if (trialNumList == 'all-at-once'):
                    data = self[r][c].getAllTrialsAsDataSet().data
                    retVar[r][c] = funReduce(getDictData(data, varList))
                else:
                    if (len(self[r][c]) == 0):
                        retVar[r, c, :] = np.nan
                    else:
                        for trialNum in trialNumList:
                            data = self[r][c][trialNum].data
                            retVar[r][c][trialNum] = funReduce(getDictData(data,
                                varList))

        return retVar
        

    @property
    def rootDir(self):
        return self._rootDir

