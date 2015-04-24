''''Classes to define abstractions on parameter exploration data.'''
from __future__ import absolute_import
from collections    import Sequence
from os.path        import exists, basename
import subprocess

import numpy as np
from ..otherpkg.log   import log_warn, log_info, getClassLogger

from ..data_storage       import DataStorage
from ..data_storage.dict  import getDictData
from .data_sets           import DictDataSet

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

    Parameters
    ----------
    fileName : str
        File path that contains the data for all the trials.
    fileMode : str
        File open mode.
    data_set_cls : class object
        A DataSet class object that will be created when accessing the data.
    '''

    def __init__(self, fileName, fileMode, data_set_cls=None):
        self._fileName = fileName
        self._fileMode = fileMode
        self._dataLoaded = False
        if data_set_cls is None:
            self._data_set_cls = DictDataSet
        else:
            self._data_set_cls = data_set_cls
        DataSpace.__init__(self, None, key='trials')

    @property
    def file_path(self):
        '''Return the full path to the file name for this trial set.'''
        return self._fileName

    @property
    def file_name_base(self):
        '''Return just the file name for this trial set.'''
        return basename(self._fileName)

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
        return self._data_set_cls(self._vals[key])

    def getAllTrialsAsDataSet(self):
        self._loadData()
        return self._data_set_cls(self._ds)

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
    '''A 2D parameter sweep space with a number of trials per job.

    .. todo::
        This class will need a massive refactoring indeed!

    Parameters
    ----------
    shape : a pair of ints, or None
        Parameter space shape (rows, columns). If ``None``, the shape will be
        retrieved automatically from the metadata file.
    rootDir : str
        Root directory for the space.
    dataPoints : a list of pairs, or None
        If this is not ``None``, then it must contain a list of coordinates in
        the data space to restrict the data manipulation to.
    fileMode : str
        File mode to open the jobs in.
    fileFormat : str
        A template for job file name formatting.
    forceWMode : bool
        Whether to force write mode.
    checkParams : bool
        Whether to check integrity of iterated parameters.
    metadata_extractor : MetaDataExtractor, or None
        Extractor for the iteration metadata, i.e. which parameters have been
        iterated, their labels, etc. If ``None``, then it will not be possible
        to extract the iteration metadata for this parameter space.
    '''
    def __init__(self, shape, rootDir, dataPoints=None, fileMode='r+',
                 fileFormat="job{0:05}_output.h5", forceWMode=False,
                 checkParams=False, metadata_extractor=None):
        if fileMode == 'w' and forceWMode == False:
            raise ValueError("'w' file open mode is not allowed. Use "
                             "'forceWMode' to override.")
        self._fileMode = fileMode
        self._rootDir = rootDir
        self._iter_file = None
        self._shape = self._determine_shape(shape)
        self._dataPoints = dataPoints
        if self._dataPoints is not None:
            self._partial = True
        else:
            self._partial = False
        self._fileFormat = fileFormat
        self._checkParams = checkParams
        self.rows = self.shape[0]
        self.cols = self.shape[1]
        self._loadTrials()
        self._aggregationDS = None
        self.saveDataFileName = 'reductions.h5'
        self._extractor = metadata_extractor

    @property
    def _meta_file(self):
        '''The iteration metadata file for this space.'''
        if self._iter_file is None:
            self._iter_file = self._open_iter_file()
        return self._iter_file

    def _open_iter_file(self):
        '''Open the iteration metadata file for this parameter space.'''
        fname = "{0}/iterparams.h5".format(self._rootDir)
        try:
            ds = DataStorage.open(fname, 'r')
        except IOError as e:
            job2DLogger.error("Could not open the metadata file. Check that "
                              "the file exists:\n\t%s", fname)
            raise e
        return ds

    @property
    def metadata(self):
        '''Return a reference to the metadata associated with this parameter
        space.
        '''
        return self._extractor

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
        if self._aggregationDS is not None:
            self._aggregationDS.close()
        # self._meta_file.close()

    def getShape(self):
        '''Return the shape of this parameter space.'''
        return self._shape

    @property
    def shape(self):
        '''Return the shape of the parameter space.'''
        return self._shape

    def _determine_shape(self, custom_shape):
        '''Determine the type of shape that should be used.

        Parameters
        ----------
        custom_shape : a tuple of ints or None
            If ``None``, will try to extract the shape from the metadata saved
            with the jobs. Otherwise will use the value in this parameters.

        Returns
        -------
        shape : a tuple of ints
            Returns the correct shape or raises RuntimeError if it cannot be
            determined.
        '''
        if custom_shape is not None:
            return custom_shape

        # Here, require that _meta_file contains the appropriate fields
        try:
            dims = self._meta_file['dimensions']
            return (dims[0], dims[1])
        except KeyError:
            raise LookupError('Could not retrieve the dimensions of this '
                              'space from metadata. You are probably using '
                              'an older version of data. In that case, the '
                              'shape of the parameter space must be specified '
                              'explicitly.')

    def getIteratedParameters(self, name_list=None):
        '''Retrieve the iterated parameters.

        Parameters
        ----------
        name_list : list of strings, or ``None``
            The list of iteration parameters. Can only be left ``None`` if the
            metadata of the space contains the iteration labels.

        .. deprecated::
            Use `~get_iterated_parameter` instead

        Returns
        -------
        parameters : list of 2D arrays
            Returns a list of the iterated parameters in the format [Rows,
            Columns]. If you need to retrieve the iteration labels, use
            :meth:`~get_iteration_labels`
        '''
        if name_list is None:
            name_list = self.get_iteration_labels()

        if len(name_list) != 2:
            raise ValueError("'name_list' parameter must contain exactly 2 "
                             "elements.")

        ret = []
        for nm in name_list:
            if nm is not None:
                ret.append(np.reshape(self._meta_file['iterParams'][nm],
                                    self._shape))
                if self._checkParams:
                    self._checkIteratedParameters(nm, ret[-1])
        return ret

    def get_iterated_parameter(self, dim):
        '''Return the iterated parameter data for dimension ``dim``.'''
        label = self.get_iteration_labels()[dim]
        if label is not None:
            return np.reshape(self._meta_file['iterParams'][label],
                              self._shape)
        else:
            return None


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

    def get_iteration_labels(self):
        '''Retrieve the iteration labels as a tuple of strings.'''
        try:
            return self._meta_file['dimension_labels']
        except KeyError:
            raise LookupError('Could not retrieve the iteration labels of '
                              'this space from metadata. You are probably '
                              'using an older version of data, in which case '
                              'you cannot use this method.')

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


class JobTrialSpace1D(JobTrialSpace2D):
    '''A 1D parameter sweep space with a number of trials per job.

    .. note::
        This parameter space behaves like the :class:`~JobTrialSpace2D`, but he
        number of rows will be forced to be strictly 1. Thus when indexing, it
        is necessary to use notation ``space[0][idx]``.

    .. note::
        Unlike :class:`~JobTrialSpace2D`, the shape here will be determined
        automatically from the metadata file.

    .. todo::
        If this class is going to be used extensively, it is necessary to write
        unit tests!

    Parameters
    ----------
    shape : None
        Must be set to ``None``. For backward compatibility only.
    rootDir : str
        Root directory for the space.
    kwargs : keyword arguments
        Keyword arguments passed on to :class:`~JobTrialSpace2D`.
    '''
    def __init__(self, shape, rootDir, **kwargs):
        if shape is not None:
            raise TypeError('shape must be None.')
        super(JobTrialSpace1D, self).__init__(None, rootDir, **kwargs)

    def _determine_shape(self, custom_shape):
        '''Determine the size of this 1D parameter space.

        .. note::
            The actual shape will be represented as a 2D shape: (1, size).
        '''
        if custom_shape is not None:
            return custom_shape

        dims = self._meta_file['dimensions']
        if len(dims) != 1:
            raise TypeError('You are trying to open a 1D parameter space, '
                            'but the actual data contains more than 1 '
                            'dimension (%dD).' % len(dims))
        return (1, dims[0])

    def _checkIteratedParameters(self, paramStr, toCheck):
        raise RuntimeError('This method cannot be called from within a 1D '
                           'dataset.')

    def get_iteration_labels(self):
        try:
            labels = self._meta_file['dimension_labels']
            return (None, labels[0])
        except KeyError:
            raise LookupError('Could not retrieve the iteration labels of '
                              'this space from metadata. You are probably '
                              'using an older version of data, in which case '
                              'you cannot use this method.')
