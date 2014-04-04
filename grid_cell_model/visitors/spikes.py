'''
Visitors that perform (raw) spikes analysis.
'''
import numpy as np

import analysis.spikes as aspikes
from .interface import DictDSVisitor
from otherpkg.log     import getClassLogger

import logging
FRLogger = getClassLogger("FiringRateVisitor", __name__)

__all__ = ['FiringRateVisitor', 'SpikeTrainXCVisitor', 'SpikeStatsVisitor']


class FiringRateVisitor(DictDSVisitor):
    '''
    Determine various firing rate statistics of a population of neurons on the
    spiking data dataset:
        * Average firing rate of all the neurons.

    '''

    def __init__(self, winLen, winDt, thetaStartT=None, tEnd=None, forceUpdate=False):
        '''
        Initialize the visitor.

        Parameters
        ----------
        winLen : float (ms)
            Length of the firing rate window as a fraction of the theta cycle
            time (0, 1>.
        winDt : float (ms)
            ``dt`` of the firing rate window.
        thetaStartT : float (ms)
            Start time of the theta signal. The center of the firing rate
            window will be in the middle of the theta signal. Therefore it is
            up to the user to ensure that the peak of the theta signal is in
            the middle. If None, extract from the data when performing analysis
        tEnd : float (ms)
            Analysis end time. If None, extract from the data
        forceUpdate : boolean, optional
            Whether to do the data analysis even if the data already exists.
        '''
        self.thetaStartT = thetaStartT
        self.tEnd        = tEnd
        self.winLen      = winLen
        self.winDt       = winDt
        self.forceUpdate = forceUpdate


    def _getSpikeTrain(self, data, monName, dimList):
        senders, times, N = DictDSVisitor._getSpikeTrain(self, data, monName,
                dimList)
        return aspikes.PopulationSpikes(N, senders, times)

    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        thetaStartT = self._checkAttrIsNone(self.thetaStartT,
                'theta_start_t', data)
        tEnd = self._checkAttrIsNone(self.tEnd, 'time', data)

        eSp = self._getSpikeTrain(data, 'spikeMon_e', ['Ne_x', 'Ne_y'])
        iSp = self._getSpikeTrain(data, 'spikeMon_i', ['Ni_x', 'Ni_y'])

        if (not self.folderExists(a, ['FR_e']) or self.forceUpdate):
            FRLogger.info("Analysing (FR_e)")
            eFR = eSp.avgFiringRate(thetaStartT, tEnd)
            a['FR_e'] = {
                    'all'             : eFR,
                    'avg'             : np.mean(eFR),
            }
        else:
            FRLogger.info("Data present (FR_e), skipping.")


        import pdb; pdb.set_trace()
        frE = a['FR_e']
        if not ('popSliding' in frE.keys() and 'popSlidingTimes' in
                frE.keys()):
            FRLogger.info("Analysing (sliding FR_e)")
            eSlidingFR, eSlidingFRt = eSp.slidingFiringRate(
                                        thetaStartT, tEnd, self.winDt,
                                        self.winLen)
            frE.update({
                    'popSliding'      : np.mean(eSlidingFR, axis=0),
                    'popSlidingTimes' : eSlidingFRt
            })
        else:
            FRLogger.info("Data present (sliding FR_e), skipping.")


        if (not self.folderExists(a, ['FR_i']) or self.forceUpdate):
            FRLogger.info("Analysing (FR_i)")
            iFR = iSp.avgFiringRate(thetaStartT, tEnd)
            a['FR_i'] = {
                    'all'             : iFR,
                    'avg'             : np.mean(iFR)
            }
        else:
            FRLogger.info("Data present (FR_i), skipping.")


        frI = a['FR_i']
        if not ('popSliding' in frI.keys() and 'popSlidingTimes' in
                frI.keys()):
            FRLogger.info("Analysing (sliding FR_i)")
            iSlidingFR, iSlidingFRt = iSp.slidingFiringRate(
                                        thetaStartT, tEnd, self.winDt,
                                        self.winLen)
            frI.update({
                    'popSliding'      : np.mean(iSlidingFR, axis=0),
                    'popSlidingTimes' : iSlidingFRt
            })
        else:
            FRLogger.info("Data present (sliding FR_i), skipping.")


        
        


##############################################################################

class SpikeTrainXCVisitor(DictDSVisitor):
    '''
    Compute spike train crosscorrelations between spikes of neuron of a
    population.
    '''

    def __init__(self, monitorName, bins, lagRange=None, neuronIdx=None,
            forceUpdate=False):
        '''
        Parameters:

        monitorName : string
            Name of the monitor in the data hierarchy.
        bins : int
            Number of bins for the cross-correlation histogram. Bin centers
        lagRange : (tStart, tEnd)
            Start/end time of the lags in histograms. This values will define
            the left and right edges of the histogram. All the values outside
            this range will be ignored.
        neuronIdx : list of ints or None
            List of neurons for which to compute the pairwise
            cross-correlation. If None, use all the neurons.
        '''
        self.allowedMonitors = ['spikeMon_e', 'spikeMon_i']
        if (not monitorName in self.allowedMonitors):
            msg = "monitorName must be one of {0}".format(allowedMonitors)
            raise ValueError(msg)

        self.monitorName = monitorName
        self.lagRange    = lagRange
        self.bins        = bins
        self.forceUpdate = forceUpdate
        self.neuronIdx   = neuronIdx

        if (self.monitorName == "spikeMon_e"):
            self.NName = "net_Ne"
            self.outputName = "XCorrelation_e"
        elif (self.monitorName == "spikeMon_i"):
            self.NName = "net_Ni"
            self.outputName = "XCorrelation_i"


    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (self.outputName in a.keys() and not self.forceUpdate):
            log_info("SpikeTrainXCorrelation", "Data present. Skipping analysis.")
            return

        spikes = simei.MonitoredSpikes(data, self.monitorName, self.NName)
        if (self.neuronIdx is None):
            idx1 = range(spikes.N)
        else:
            idx1 = self.neuronIdx
        correlations, bin_centers, bin_edges = spikes.spikeTrainXCorrelation(idx1, None,
                self.lagRange, self.bins)

        a[self.outputName] = dict(
                neuronIdx    = idx1,
                correlations = correlations,
                bin_edges    = bin_edges,
                bin_centers  = bin_centers)


class SpikeStatsVisitor(DictDSVisitor):
    def __init__(self, monitorName, forceUpdate=False):
        '''
        Parameters:

        monitorName : string
            Name of the monitor in the data hierarchy.
        '''
        self.allowedMonitors = ['spikeMon_e', 'spikeMon_i']
        if (not monitorName in self.allowedMonitors):
            msg = "monitorName must be one of {0}".format(allowedMonitors)
            raise ValueError(msg)

        self.monitorName = monitorName
        self.forceUpdate = forceUpdate

        if (self.monitorName == "spikeMon_e"):
            self.NName = "net_Ne"
            self.outputName = "CV_e"
        elif (self.monitorName == "spikeMon_i"):
            self.NName = "net_Ni"
            self.outputName = "CV_i"


    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (self.outputName in a.keys() and not self.forceUpdate):
            log_info("SpikeStatsVisitor", "Data present. Skipping analysis.")
            return

        spikes = simei.MonitoredSpikes(data, self.monitorName, self.NName)
        a[self.outputName] = np.array(spikes.ISICV())
