'''
Grid field plotting visitors.
'''
import os
import errno

import numpy as np
from matplotlib.pyplot import figure, plot, pcolormesh, subplot2grid, savefig,\
        colorbar, axis, xlabel, ylabel

from analysis.spikes import PopulationSpikes
from analysis.grid_cells import SNSpatialRate2D, SNAutoCorr, cellGridnessScore
from plotting.bumps  import torusFiringRate
from plotting.grids  import plotSpikes2D
from otherpkg.log    import log_warn, log_info

from .. import interface


__all__ = ['GridPlotVisitor']



class GridPlotVisitor(interface.DictDSVisitor):
    '''
    Plot the following figures for the data set:
        * Spike map, showing the trajectory of a simulated animal and spike
          positions
        * Smoothed rate map
        * Autocorrelation of the rate map
    '''
    class PlotOptions(object):
        def __init__(self):
            self.bump        = True
            self.spikes      = True
            self.rateMap     = True
            self.fft         = False
            self.sn_ac       = True
            self.gridness_ac = True

        def setAll(self, val):
            self.bump        = val
            self.spikes      = val
            self.rateMap     = val
            self.fft         = val
            self.sn_ac       = val
            self.gridness_ac = val

    class BumpOnly(PlotOptions):
        def __init__(self):
            GridPlotVisitor.PlotOptions.__init__(self)
            self.setAll(False)
            self.bump = True



    def __init__(self, rootDir, spikeType='E', neuronNum=0, arenaDiam=180.0,
            smoothingSigma=3.0, bumpTStart=None, bumpTEnd=None,
            minGridnessT=0.0, plotOptions=PlotOptions(), forceUpdate=False):
        '''
        Parameters
        ----------
        rootDir : string
            Root output directory where to store the output figures.
        spikeType : string
            Type of the analysis, can be either 'E' - to plot grid field data
            for E cells, or 'I' to plot them for the I cells.
        neuronNum : int, optional
            Neuron index.
        arenaDiam : float (cm)
            Arena diameter
        smoothingSigma : float (cm), optional
            Gaussian smoothing kernel width
        plotBump : bool, optional
            Whether to plot bump (needs spiking activity of the whole
            population)
        bumpTStart : float
            Start time for bump (firing rate) plots (ms)
        bumpTEnd   : float
            End time for bump plots (ms)
        minGridnessT : float
            Minimal time (determined by the time of last spike of an E cell) to
            consider the simulation for gridness score (if less than this time,
            gridness score will be NaN).
        plotOptions : PlotOptions
            A set of flags that specify what should be plotted. TODO: this is
            broken now, since there are dependencies.
        forceUpdate : bool
            Whether to force data analysis and saving of the results even when
            they already exist.
        '''
        self.rootDir        = rootDir
        self.neuronNum      = neuronNum
        self.outputDir      = 'grids'
        self.arenaDiam      = arenaDiam
        self.smoothingSigma = smoothingSigma
        self.bumpTStart     = bumpTStart
        self.bumpTEnd       = bumpTEnd
        self.po             = plotOptions
        self.minGridnessT   = minGridnessT
        self.forceUpdate    = forceUpdate
        self.setSpikeType(spikeType)



    def _checkSpikeType(self, t):
        if (t == 'E' or t == 'I'):
            return True
        msg = "spikeType must be 'E' or 'I'. Got '{0}'".format(spikeType)
        raise ValueError(msg)

    def setSpikeType(self, t):
        self._checkSpikeType(t)
        self._spikeType = t


    def createOutputDirs(self):
        # Create output directory(ies)
        try:
            os.makedirs('{0}/{1}'.format(self.rootDir, self.outputDir))
        except OSError as e:
            if (e.errno == errno.EEXIST):
                log_warn("GridPlotVisitor", 'Output directory already ' +
                        'exists. This might overwrite files.')
            else:
                raise e


    
    def _shiftSpikeTimes(self, s, startT):
        '''
        Shift the times of the spikes so that the time starts at startT. All
        the spikes with negative time will be removed
        '''
        s -= startT
        return np.delete(s, np.nonzero(s < 0)[0])


    def visitDictDataSet(self, ds, **kw):
        data = ds.data

        if ('analysis' in data.keys() and not self.forceUpdate):
            log_info("GridPlotVisitor", "Data present. Skipping analysis.")
            return

        simT = self.getOption(data, 'time') # ms
        jobNum = self.getOption(data, 'job_num')
        if ('trialNum' in kw.keys()):
            trialNum = kw['trialNum']
        else:
            trialNum = 0
        self.createOutputDirs()
        fileNameTemplate = "{0}/{1}/job{2:05}_trial{3:03}".format(self.rootDir,
                self.outputDir, jobNum, trialNum)

        pos_x         = self.getNetParam(data, 'rat_pos_x')
        pos_y         = self.getNetParam(data, 'rat_pos_y')
        rat_dt        = self.getNetParam(data, 'rat_dt')
        velocityStart = self.getOption(data, 'theta_start_t')
        if self._spikeType == 'E':
            monName = 'spikeMon_e'
        if self._spikeType == 'I':
            monName = 'spikeMon_i'

        gridSep = self.getOption(data, 'gridSep')
        corr_cutRmin = gridSep / 2

        spikes = DictDSVisitor._getNeuronSpikeTrain(data, monName, self.neuronNum)
        spikes = self._shiftSpikeTimes(spikes, velocityStart)

        out = {}
        
        if (self.po.bump):
            figure()
            if (self.bumpTStart is None):
                bumpTStart = velocityStart
            if (self.bumpTEnd is None):
                bumpTEnd = bumpTStart + 1e3
            senders, times, N = DictDSVisitor._getSpikeTrain(self, data, monName,
                    ['Ne_x', 'Ne_y'])
            sp = PopulationSpikes(N, senders, times)
            Fe = sp.avgFiringRate(bumpTStart, bumpTEnd)
            Ne_x = self.getNetParam(data, 'Ne_x')
            Ne_y = self.getNetParam(data, 'Ne_y')
            bump_e = np.reshape(Fe, (Ne_y, Ne_x))
            torusFiringRate(
                    rateMap  = bump_e,
                    labelx   = '',
                    labely   = 'Neuron #',
                    titleStr = 'E firing rate')
            savefig(fileNameTemplate + '_bump_' + self._spikeType + '.png')
            out['bump_e'] = bump_e


        if (self.po.spikes):
            figure()
            plotSpikes2D(spikes, pos_x, pos_y, rat_dt)
            savefig(fileNameTemplate + '_spikePlot_' + self._spikeType + '.png')
            out['spikes_e'] = spikes
            out['rat_pos_x'] = pos_x
            out['rat_pos_y'] = pos_y
            out['rat_dt']    = rat_dt

        if (self.po.rateMap):
            figure()
            rateMap, xedges, yedges = SNSpatialRate2D(spikes, pos_x, pos_y, rat_dt,
                    self.arenaDiam, self.smoothingSigma)
            rateMap *= 1e3 # should be Hz
            X, Y = np.meshgrid(xedges, yedges)
            pcolormesh(X, Y, rateMap)
            colorbar()
            axis('equal')
            axis('off')
            savefig('{0}_rateMap_{1}.png'.format(fileNameTemplate,
                self._spikeType))
            out['rateMap_e'] = rateMap
            out['rateMap_e_X'] = X
            out['rateMap_e_Y'] = Y

        
        if (self.po.fft):
            figure()
            FT_size = 256
            Fs = 1.0/(self.smoothingSigma/100.0) # units: 1/m
            rateMap_pad = np.ndarray((FT_size, FT_size))
            rateMap_pad[:, :] = 0
            rateMap_pad[0:rateMap.shape[0], 0:rateMap.shape[0]] = rateMap - np.mean(rateMap.flatten())
            FT = fft2(rateMap_pad)
            fxy = np.linspace(-1.0, 1.0, FT_size)
            fxy_igor = Fs/2.0*np.linspace(-1.0, 1.0, FT_size+1)
            FX, FY = np.meshgrid(fxy, fxy)
            FX *= Fs/2.0
            FY *= Fs/2.0
            PSD_centered = np.abs(np.fft.fftshift(FT))**2
            pcolormesh(FX, FY, PSD_centered)
            #axis('equal')
            xlim([-10, 10])
            ylim([-10, 10])
            savefig('{0}_fft2_{1}.png'.format(fileNameTemplate, self._spikeType))
            out['FFTX'] = FX
            out['FFTY'] = FY
            out['FFT']  = PSD_centered


        if (self.po.sn_ac):
            figure()
            corr, xedges_corr, yedges_corr = SNAutoCorr(rateMap, self.arenaDiam,
                self.smoothingSigma)
            X, Y = np.meshgrid(xedges_corr, yedges_corr)
            pcolormesh(X, Y, corr)
            axis('equal')
            axis('off')
            savefig('{0}_rateCorr_{1}.png'.format(fileNameTemplate,
                self._spikeType))
            out['corr_X'] = X
            out['corr_Y'] = Y
            out['corr']   = corr


        if (self.po.gridness_ac):
            figure()
            G, crossCorr, angles = cellGridnessScore(rateMap, self.arenaDiam,
                self.smoothingSigma, corr_cutRmin)
            plot(angles, crossCorr)
            xlabel('Angle (deg.)')
            ylabel('Corr. coefficient')
            savefig('{0}_gridnessCorr_{1}.png'.format(fileNameTemplate,
                self._spikeType))
            # Gridness score valid only when T >= minGridnessT
            lastSpikeT = data['spikeMon_e']['events']['times'][-1]
            if (lastSpikeT >= self.minGridnessT):
                out['gridnessScore']  = G
            else:
                log_warn('GridPlotVisitor', 'Simulation too short, G <- NaN')
                out['gridnessScore']  = np.nan
            out['gridnessCorr']   = crossCorr
            out['gridnessAngles'] = angles

        data['analysis'] = out



