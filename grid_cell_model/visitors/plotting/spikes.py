'''
Visitors that perform plotting of spikes.
'''
import numpy as np


from otherpkg.log import getClassLogger
#import analysis.spikes as aspikes
from plotting.signal import signalPlot

from .. import interface
from .. import spikes

import logging
FRPlotLogger = getClassLogger("FiringRatePlotter", __name__)


class FiringRatePlotter(interface.DictDSVisitor):
    '''Plot population firing rates for the specified duration'''

    def __init__(self, rootDir=None, readme=''):
        '''Initialize the visitor

        Parameters
        ----------
        rootDir : str
            Root output directory where to save the image. It will be a
            subdirectory of where the data is located.
        '''
        super(FiringRatePlotter, self).__init__()
        self.rootDir = rootDir



    def visitDictDataSet(self, ds, **kw):
        if 'fileName' not in kw.keys():
            msg = 'Did not receive the fileName as a keyword argument.'
            velGainLogger.warn(msg)
            return
        r = kw.pop('r', '')
        c = kw.pop('c', '')
        trialNum = kw.pop('trialNum', '')

        data = ds.data
        a = data['analysis']
        fig = plt.Figure()
        
        # E firing rate
        FR_e  = a['FR_e/popSliding']
        FRt_e = a['FR_e/popSligindTimes']
        
        axE = fig.add_subplot(211)
        signalPlot(FRt_e, FR_e, axE) 


        # I firing rate
        FR_i  = a['FR_i/popSliding']
        FRt_i = a['FR_i/popSligindTimes']
        
        axI = fig.add_subplot(212)
        signalPlot(FRt_i, FR_i, axI) 

        figPath = self.getFigPath(kw['fileName'], self.rootDir, r, c, trialNum)
        FRPlotLogger.info("Saving figure to '{0}'".format(figPath))
        fig.savefig(figPath)


################################################################################
#class ISIPlotVisitor(DictDSVisitor):
#    def __init__(self, rootDir, spikeType, nCols, nRows, **kw):
#        '''
#        Parameters
#        ----------
#        rootDir : string
#            Root output directory where to store the output figures.
#        spikeType : string
#            Type of the analysis, can be either 'E' - to plot grid field data
#            for E cells, or 'I' to plot them for the I cells.
#        nCols : int
#            Number of columns in the grid plot
#        nRows : int
#            Number of rows in the grid plot
#        '''
#        self.rootDir   = rootDir
#        self.outputDir = 'ISI_statistics'
#        self.setSpikeType(spikeType)
#        self.nCols     = nCols
#        self.nRows     = nRows
#
#        self.ISINWindows = kw.pop('ISINWindows', 0)
#
#        self.hist_kw   = kw
#
#
#    def _checkSpikeType(self, t):
#        if (t == 'E' or t == 'I'):
#            return True
#        msg = "spikeType must be 'E' or 'I'. Got '{0}'".format(spikeType)
#        raise ValueError(msg)
#
#    def setSpikeType(self, t):
#        self._checkSpikeType(t)
#        self._spikeType = t
#
#
#    def createOutputDirs(self):
#        # Create output directory(ies)
#        try:
#            os.makedirs('{0}/{1}'.format(self.rootDir, self.outputDir))
#        except OSError as e:
#            if (e.errno == errno.EEXIST):
#                log_warn("GridPlotVisitor", 'Output directory already ' +
#                        'exists. This might overwrite files.')
#            else:
#                raise e
#
#
#    def visitDictDataSet(self, ds, **kw):
#        data = ds.data
#
#        simT = self.getOption(data, 'time') # ms
#        thetaT = 1e3 / self.getOption(data, 'theta_freq')
#        jobNum = self.getOption(data, 'job_num')
#        if ('trialNum' in kw.keys()):
#            trialNum = kw['trialNum']
#        else:
#            trialNum = 0
#        self.createOutputDirs()
#        fileNameTemplate = "{0}/{1}/job{2:05}_trial{3:03}_{4}".format(self.rootDir,
#                self.outputDir, jobNum, trialNum, self._spikeType)
#
#        if self._spikeType == 'E':
#            monName = 'spikeMon_e'
#            NName   = 'net_Ne'
#        if self._spikeType == 'I':
#            monName = 'spikeMon_i'
#            NName   = 'net_Ni'
#
#         
#        # Pick the most-spiking neurons
#        spikes = MonitoredSpikes(data, monName, NName)
#        rate = spikes.avgFiringRate(0, simT)
#        maxRateIdx = np.argsort(rate)
#        ISIs = spikes.ISI(maxRateIdx[0:self.nRows*self.nCols])
#        
#        ## ISI histogram plots
#        #fig = plt.figure(figsize=(11.69, 6.57))
#        #gs = plt.GridSpec(self.nRows, self.nCols)
#        #it = 0
#        #for r in xrange(self.nRows):
#        #    for c in xrange(self.nCols):
#        #        if (it >= spikes.N):
#        #            break
#
#        #        ax = plt.subplot(gs[r, c]) 
#        #        ax.hist(ISIs[it], **self.hist_kw)
#        #        if (r == self.nRows - 1):
#        #            ax.set_xlabel('ISI (ms)')
#        #        ax.xaxis.set_major_locator(ti.LinearLocator(2))
#        #        ax.set_yticks([])
#        #        it += 1
#
#        #gs.tight_layout(fig)
#        #fname = '{0}_isi_histograms.pdf'.format(fileNameTemplate)
#        #fig.savefig(fname)
#
#        # CV plots, as a function of window length
#        if (self.ISINWindows != 0):
#            winLens = np.arange(1, self.ISINWindows+1)*thetaT
#            CVs = spikes.ISICV(n=maxRateIdx[0:self.nRows*self.nCols],
#                    winLen=winLens)
#            fig = plt.figure(figsize=(11.69, 6.57))
#            gs = plt.GridSpec(self.nRows, self.nCols)
#            it = 0
#            for r in xrange(self.nRows):
#                for c in xrange(self.nCols):
#                    if (it >= spikes.N):
#                        break
#
#                    ax = plt.subplot(gs[r, c]) 
#                    ax.plot(winLens, CVs[it])
#                    if (r == self.nRows - 1):
#                        ax.set_xlabel('Window (ms)')
#                    ax.xaxis.set_major_locator(ti.LinearLocator(2))
#                    ax.yaxis.set_major_locator(ti.MaxNLocator(3))
#                    it += 1
#
#            gs.tight_layout(fig)
#            fname = '{0}_isi_CV_window.pdf'.format(fileNameTemplate)
#            fig.savefig(fname)
#
#        
