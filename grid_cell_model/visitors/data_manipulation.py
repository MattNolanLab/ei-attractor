'''
Data manipulation visitors. Visitors that do some specific/generic
manipulations with the data they receive.
'''
import logging
from interface        import DictDSVisitor 
from otherpkg.log import getClassLogger

velPruneLogger = getClassLogger('VelocityPruningVisitor', __name__)

class VelocityPruningVisitor(DictDSVisitor):
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
        self.what = what


    def visitDictDataSet(self, ds, **kw):
        if ds.data is None:
            velPruneLogger.info('No data set, skipping everything.')
            return
        data = ds.data
        trials = data['trials']

        for trialNum in xrange(len(trials)):
            velPruneLogger.info("Trial no. %d/%d", trialNum, len(trials)-1)
            IvelVec = trials[trialNum]['IvelVec']
            for IvelIdx in xrange(len(IvelVec)):
                iData = trials[trialNum]['IvelData'][IvelIdx]
                velPruneLogger.info("Processing vel. index: %d/%d; %.2f pA",
                        IvelIdx, len(IvelVec) - 1, IvelVec[IvelIdx])

                if self.what in iData.keys():
                    del iData[self.what]
                    velPruneLogger.info("'%s' deleted.", self.what)
                else:
                    velPruneLogger.info("'%s' does not exist. No delete",
                            self.what)

