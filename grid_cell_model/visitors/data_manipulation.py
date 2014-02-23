'''
Data manipulation visitors. Visitors that do some specific/generic
manipulations with the data they receive.
'''
from abc import abstractmethod
import logging

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
    Convert specified data in the data set.
    '''

    def __init__(self, what, targetType):
        self.what = what
        self.targetType = targetType

    def manipulateData(self, data, **kw):
        if self.what not in data.keys():
            convLogger.warn("'%s' not in present in data. Quitting.",
                    self.what)
            return

        try:
            tmp = self.targetType(data[self.what])
        except e:
            convLogger.warn('Conversion failed. The original data is intact')
            return

        data[self.what] = tmp




