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
from collections    import Sequence
from os.path        import exists
import subprocess

import numpy as np
from otherpkg.log   import log_warn, log_info, getClassLogger

from data_storage       import DataStorage
from data_storage.dict  import getDictData
from data_sets          import DictDataSet

__all__ = []


job2DLogger = getClassLogger('JobTrialSpace2D', __name__)

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



trialSetLogger = getClassLogger('TrialSet', __name__)
class TrialSet(DataSpace):
    '''
    A 1D DataSpace that contains a list of DataSet objects.
    '''

    def __init__(self, fileName, fileMode):
        self._fileName = fileName
        self._fileMode = fileMode
        self._dataLoaded = False
        DataSpace.__init__(self, None, key='trials')


    def _loadData(self):
        if (self._dataLoaded):
            return
        try:
            trialSetLogger.debug("Opening " + self._fileName)
            self._ds = DataStorage.open(self._fileName, self._fileMode)
            DataSpace.__init__(self, self._ds['trials'], key='trials')
        except (IOError, KeyError) as e:
            self._ds = None
            msg =  "Could not open file {0}. Creating an empty DataSet instead."
            trialSetLogger.warn(msg.format(self._fileName))
            trialSetLogger.warn("Error message: {0}".format(str(e)))
            DataSpace.__init__(self, [], key='trials')
        self._dataLoaded = True


    def __len__(self):
        self._loadData()
        return DataSpace.__len__(self)


    def __del__(self):
        if (self._dataLoaded):
            if (self._ds is not None):
                trialSetLogger.debug("Closing: %s", self._fileName)
                self._ds.close()

    def __getitem__(self, key):
        self._loadData()
        return DictDataSet(self._vals[key])

    def getAllTrialsAsDataSet(self):
        self._loadData()
        return DictDataSet(self._ds)

    def visit(self, visitor, trialList=None, **kw):
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
            self.getAllTrialsAsDataSet().visit(visitor,
                                               fileName=self._fileName, **kw)
        else:
            if (trialList is None):
                trialList = xrange(len(self))

            for trialIdx in trialList:
                trial = self[trialIdx]
                trial.visit(visitor, fileName=self._fileName,
                            trialNum=trialIdx, **kw)


class DummyTrialSet(DataSpace):
    '''An empty DataSpace that does nothing.'''
    def __init__(self):
        DataSpace.__init__(self, [], key='trials')

    def __getitem__(self, key):
        raise ValueError("A DummyTrialSet cannot be indexed!")

    def visit(self, visitor, trialList=None, **kw):
        pass



class JobTrialSpace2D(DataSpace):
    def __init__(self, shape, rootDir, dataPoints=None, fileMode='r+',
            fileFormat="job{0:05}_output.h5", forceWMode=False,
            checkParams=False):
        '''
        Initialize this parameter space object.
        TODO: this class will need a massive refactoring indeed!

        Parameters
        ----------
        TODO
        '''
        if (fileMode == 'w' and forceWMode == False):
            raise ValueError("'w' file open mode is not allowed. Use " +
                    "'forceWMode' to override.")
        self._fileMode = fileMode
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
        self._loadTrials()
        self._aggregationDS = None
        self.saveDataFileName = 'reductions.h5'


    def repackItem(self, r, c):
        '''
        Repack the underlying item at ``r`` and ``c`` position.

        .. warning::

            This applies only to HDF5, use with care.
        '''

        def _cleanUp(r, c, tmpFile):
            job2DLogger.info('Cleaning up temporary.')
            subprocess.call(['rm', '-f', tmpFile])
            job2DLogger.info("Reloading the new item")
            self._vals[r]._vals[c] = self._loadItem(r, c)

        job2DLogger.info("Repacking item at %d, %d", r, c)
        currentFile = self._getFilename(r, c)
        tmpFile = currentFile + ".repacked"


        job2DLogger.info("Removing current item")
        #import pdb; pdb.set_trace()
        self._vals[r]._vals[c] = None

        # Try to repack to temporary file
        job2DLogger.info("Repacking to temporary file...")
        repackArgs = ['h5repack', currentFile, tmpFile]
        if subprocess.call(repackArgs) != 0:
            job2DLogger.warn('Repacking to %s failed.', tmpFile)
            _cleanUp(r, c, tmpFile)
            return

        # Move temporary file to original
        job2DLogger.info("Replacing original file")
        moveArgs = ['mv', tmpFile, currentFile]
        if subprocess.call(moveArgs) != 0:
            job2DLogger.warn('Replacing original %s failed.', tmpFile)
            _cleanUp(r, c, tmpFile)
            return

        # Reload the repacked file
        _cleanUp(r, c, tmpFile)


    def repackAllItems(self):
        for row in xrange(self.shape[0]):
            for col in xrange(self.shape[1]):
                self.repackItems(row, col)


    def _getFilename(self, row, col):
        it = row * self.shape[1] + col
        fileName = self._rootDir + '/' + self._fileFormat.format(it)
        #job2DLogger.debug('row: %d, col: %d, filename: %s', row, col, fileName)
        return fileName
            

    def _loadItem(self, row, col):
        # Here either we have a full data space, or a particular row
        # and column are present in the list of selected data points.
        # Otherwise, append an empty list, i.e. this will do nothing
        if (self._dataPoints is None or
                (row, col) in self._dataPoints):
            fileName = self._getFilename(row, col)
            return TrialSet(fileName, self._fileMode)
        else:
            return DummyTrialSet()


    def _loadTrials(self):
        '''
        Initialize the data space fully. Opening all the data files as
        specified in the constructor.
        '''
        rowData = []
        for row in xrange(self.rows):
            colData = []
            for col in xrange(self.cols):
                colData.append(self._loadItem(row, col))
            rowData.append(DataSpace(colData))
        DataSpace.__init__(self, rowData)


    def __len__(self):
        return self._shape[0] 

    def __del__(self):
        if (self._aggregationDS is not None):
            self._aggregationDS.close()

    def getShape(self):
        return self._shape

    @property
    def shape(self):
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
                self[r][c].visit(visitor, trialList=trialList, r=r, c=c)


    def _createAggregateOutput(self, trialNumList, output_dtype):
        '''
        Create the correct data type given the parameters. The output can be
        either a numpy.ndarray or a multidimensional list.

        Parameters
        ----------
        see self.aggregateData
        '''
        shape = self.getShape()
        rows = shape[0]
        cols = shape[1]
        retVar = None
        if (output_dtype == 'array'):
            if (trialNumList == 'all-at-once'):
                retVar = np.ndarray((rows, cols))
            else:
                nTrials = len(trialNumList)
                retVar = np.ndarray((rows, cols, nTrials))
        elif (output_dtype == 'list'):
            if (trialNumList == 'all-at-once'):
                retVar = [[np.nan for x in xrange(cols)] for x in xrange(rows)]
            else:
                nTrials = len(trialNumList)
                retVar = [[[np.nan for t in xrange(nTrials)] for c in
                    xrange(cols)] for r in xrange(rows)]
        return retVar


    def _aggregateItem(self, retVar, r, c, trialNumList, varList, funReduce):
        if (trialNumList == 'all-at-once'):
            data = self[r][c].getAllTrialsAsDataSet().data
            if (data is None):
                retVar[r][c] = np.nan
            else:
                try:
                    retVar[r][c] = funReduce(getDictData(data, varList))
                except (IOError, KeyError) as e:
                    self._reductionFailureMsg(e, r, c)
                    retVar[r][c] = np.nan
        else:
            for trialNum in trialNumList:
                if (len(self[r][c]) == 0):
                    retVar[r][c][trialNum] = np.nan
                    continue

                try:
                    data = self[r][c][trialNum].data
                    retVar[r][c][trialNum] = funReduce(getDictData(data,
                        varList))
                except (IOError, KeyError) as e:
                    self._reductionFailureMsg(e, r, c)
                    retVar[r][c][trialNum] = np.nan


    def _reductionFailureMsg(self, e, r, c):
        msg = 'Reduction step failed at (r, c) == ({0}, '+\
            '{1}). Setting value as NaN.'
        log_warn('JobTrialSpace2D', msg.format(r, c))
        log_warn("JobTrialSpace2D", "Error message: {0}".format(str(e)))


    def _getAggregationDS(self):
        nm = '{0}/{1}'.format(self._rootDir, self.saveDataFileName)
        try:
            self._aggregationDS = DataStorage.open(nm, 'a')
            return self._aggregationDS
        except IOError as e:
            return None


    def getReduction(self, path):
        inData = self._getAggregationDS()
        if isinstance(path, str):
            return inData[path]
        else:
            return inData.getItemChained(path)


    def aggregateData(self, varList, trialNumList, funReduce=None,
            output_dtype='array', loadData=True, saveData=False,
            saveDataFileName='reductions.h5'):
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
        output_dtype : string
            Output array data type. If 'array' the output will be of type
            numpy.ndarray; if 'list', the output will be list.
        loadData : bool, optional
            If True, try to load data from the reductions data file (defined by
            self.saveDataFileName). If the data cannot be loaded, do the
            reductions itself. Setting this parameter to True imlies 'saveData'
            will be False.
        saveData : bool, optional
            Whether to save data to a file, under the rootDir directory. The
            data will be saved at the top level of a dictionary data set, with
            a key taken from the last item of varList
        output : A 3D numpy array if trialNumList is a list, or a 2D array
                 otherwise
            All the aggregated data
        '''
        if (self._partial):
            # Cannot do aggregation on a restricted data set
            raise NotImplementedError("Data aggregation on a partial data " +
                    "space has not been implemented yet.")

        # Try to load data
        if (loadData):
            try:
                msg = 'Loading aggregated data from file: {0}, vars: {1}'
                log_info('JobTrialSpace2D', msg.format(self.saveDataFileName,
                    varList))
                inData = self._getAggregationDS()
                if (inData is not None):
                    return inData.getItemChained(varList)
                else:
                    io_err = 'Could not open file: {0}. Performing the reduction.'
                    log_info('JobTrialSpace2D', io_err.format(self.saveDataFileName))
            except KeyError as e:
                key_err = 'Could not load var: {0}. Performing the reduction.'
                log_info('JobTrialSpace2D', key_err.format(varList[-1]))


        # Data not loaded, do the reduction
        shape = self.getShape()
        rows = shape[0]
        cols = shape[1]
        retVar = self._createAggregateOutput(trialNumList, output_dtype)

        if (funReduce is None):
            funReduce = lambda x: x

        for r in xrange(rows):
            for c in xrange(cols):
                self._aggregateItem(retVar, r, c, trialNumList, varList, funReduce)


        if (saveData):
            msg = 'Saving aggregated data into file: {0}, vars: {1}'
            log_info('JobTrialSpace2D', msg.format(self.saveDataFileName,
                varList))
            outData = self._getAggregationDS()
            if (outData is not None):
                outData.setItemChained(varList, retVar)
            else:
                io_err = 'Could not open file: {0}. Not saving the reduced data!'
                log_warn('JobTrialSpace2D', io_err.format(self.saveDataFileName))

        return retVar
        

    @property
    def rootDir(self):
        return self._rootDir

