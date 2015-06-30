'''Data manipulation visitors.

.. currentmodule:: grid_cell_model.visitors.data_manipulation

Visitors that do some specific/generic manipulations with the data they
receive.

Classes
-------

.. autosummary::

    VelocityDataVisitor
    VelocityPruningVisitor
    VelocityConversionVisitor
'''
import collections
from abc import abstractmethod
import logging
import numpy as np

from interface        import DictDSVisitor
from otherpkg.log import getClassLogger

dataLogger  = getClassLogger('VelocityDataVisitor', __name__)
pruneLogger = getClassLogger('VelocityPruningVisitor', __name__)
convLogger  = getClassLogger('VelocityConversionVisitor', __name__)

class VelocityDataVisitor(DictDSVisitor):

    def visitDictDataSet(self, ds, **kw):
        if ds.data is None:
            dataLogger.info('No data set, skipping everything.')
            return
        data = ds.data
        trials = data['trials']

        for trialNum in xrange(len(trials)):
            dataLogger.info("Trial no. %d/%d", trialNum, len(trials)-1)
            IvelVec = trials[trialNum]['IvelVec']
            for IvelIdx in xrange(len(IvelVec)):
                iData = trials[trialNum]['IvelData'][IvelIdx]
                dataLogger.info("Processing vel. index: %d/%d; %.2f pA",
                        IvelIdx, len(IvelVec) - 1, IvelVec[IvelIdx])
                self.manipulateData(iData)

    @abstractmethod
    def manipulateData(self, data, **kw):
        raise NotImplementedError()


class VelocityPruningVisitor(VelocityDataVisitor):
    '''
    Delete certain data on each trial, for each IvelData iteration.
    '''
    def __init__(self, what):
        '''
        **Parameters**

        what : str
            Specifies what data should be deleted. The deletion path will be
            trials/<trialNum>/IvelData/<IvelIdx>/<what>, for each trialNum and
            IvelIdx.
        '''
        self.what = what


    def manipulateData(self, data):
        if self.what in data.keys():
            del data[self.what]
            pruneLogger.info("'%s' deleted.", self.what)
        else:
            pruneLogger.info("'%s' does not exist. No delete",
                    self.what)


class VelocityConversionVisitor(VelocityDataVisitor):
    '''
    Convert specified data in the data set. Apply recursively if the
    destination is an instance of Mapping.
    '''

    def __init__(self, what, targetType):
        self.what = what
        self.targetType = targetType

        self.listPath = self.what.split('/')
        if '' in self.listPath:
            raise ValueError("%s cannot contain double '/' character!" %
                    self.what)

    def manipulateData(self, data, **kw):
        try:
            root = data.get_item_chained(self.listPath[:-1])
        except KeyError as e:
            convLogger.warn("Path to '%s' does not exist. Quitting.",
                    self.listPath[-1])
            return
        self._applyConversion(root, self.listPath[-1])

    def _applyConversion(self, dataDict, key):
        item = dataDict[key]
        if isinstance(item, collections.Mapping):
            for childKey in item.keys():
                self._applyConversion(item, childKey)
        else:
            targetCls = getattr(__builtins__, self.targetType,
                    getattr(np, self.targetType))
            try:
                newItem = targetCls(item)
            except Exception as e:
                convLogger.warn('Conversion failed. The original data is intact')
                return
            dataDict[key] = newItem


