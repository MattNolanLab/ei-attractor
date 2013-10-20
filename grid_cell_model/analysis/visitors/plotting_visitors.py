#
#   plotting_visitors.py
#
#   Visitors that primarily plot the data.
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
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, pcolormesh, subplot2grid, savefig,\
        colorbar, axis, xlabel, ylabel
import os
import errno

from interface       import DictDSVisitor
from data_storage.sim_models.ei import extractStateVariable, sumAllVariables
from plotting.signal import signalPlot
from analysis.spikes import PopulationSpikes
from analysis.grid_cells import extractSpikePositions2D, SNSpatialRate2D, \
        SNAutoCorr, cellGridnessScore
from plotting.bumps  import torusFiringRate
from plotting.grids  import plotSpikes2D
from otherpkg.log    import log_warn, log_info

__all__ = ['DetailedPlotVisitor', 'GridPlotVisitor']


class DetailedPlotVisitor(DictDSVisitor):
    '''
    Make a plot of E and I membrane potentials, currents and torus firing
    rates. Annotate the plots with the data set parameters used in the noise
    parameter sweeps.

    Save the plot to a file.
    '''

    def __init__(self, rootDir, plotT, bumpTStart=None, bumpTEnd=None):
        '''
        Initialize the visitor.

        Parameters
        ----------
        rootDir : str
            Path to the output directory
        plotT : float
            Amount of time at the end of the simulation to plot (ms)
        bumpTStart : float
            Start time for bump (firing rate) plots (ms)
        bumpTEnd   : float
            End time for bump plots (ms)
        '''
        self.rootDir    = rootDir
        self.plotT      = plotT
        self.bumpTStart = bumpTStart
        self.bumpTEnd   = bumpTEnd

        self.nIdxMiddle = 0
        self.nIdxEdge = 1

    def _getSpikeTrain(self, data, monName, dimList):
        senders, times, N = DictDSVisitor._getSpikeTrain(self, data, monName,
                dimList)
        return PopulationSpikes(N, senders, times)


    def _sliceSignal(self, t, sig, tStart, tEnd):
        idx = np.logical_and(t >= tStart, t <= tEnd)
        return t[idx], sig[idx], idx

    def _extractStateVars(self, mon, varName, plotTStart, plotTEnd):
        '''
        Extract state variables from a pair of monitors. One in the centre, the
        other one at the edge of the neural sheet.
        '''
        t, dt = extractStateVariable(mon, self.nIdxMiddle, 'times')
        sig, dt = sumAllVariables(mon, self.nIdxMiddle, varName)
        t, sigMiddle, idx = self._sliceSignal(t, sig, plotTStart, plotTEnd)
        sig, dt = sumAllVariables(mon, self.nIdxEdge, varName)
        sigEdge = sig[idx]
        return t, sigMiddle, sigEdge

    def _dictData(self, d, keyList):
        ret = d
        for key in keyList:
            ret = ret[key]
        return ret

    
    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        a = data['analysis']
        fig = figure(figsize=(14, 8))
        mon_e = data['stateMon_e']
        mon_i = data['stateMon_i']
        simT = self.getOption(data, 'time') # ms
        plotTStart = simT - self.plotT
        plotTEnd   = simT
        x_lim = [plotTStart, plotTEnd]
        if ('trialNum' in kw.keys()):
            trialNum = kw['trialNum']
        else:
            trialNum = 0

        # Create output directory(ies)
        try:
            os.makedirs(self.rootDir)
        except OSError as e:
            if (e.errno == errno.EEXIST):
                log_warn("DetailedPlotVisitor", 'Output directory already ' +
                        'exists. This might overwrite files.')
            else:
                raise e
        
        # Plot E Vm
        rows = 4
        cols = 8
        ax_Vm = subplot2grid((rows, cols), (1, 0), colspan=3)
        t, VmMiddle, VmEdge = self._extractStateVars(mon_e, ['V_m'], plotTStart,
                plotTEnd)
        signalPlot(t, VmMiddle, ax_Vm, labelx="", labely = 'E cell $V_m$',
                nThetaTicks=5)
        signalPlot(t, VmEdge, ax_Vm, labelx="", labely = 'E cell $V_m$',
                nThetaTicks=5)
        plt.xlim(x_lim)

        # Plot I Vm
        ax_Vm = subplot2grid((rows, cols), (2, 0), colspan=3)
        t, VmMiddle, VmEdge = self._extractStateVars(mon_i, ['V_m'], plotTStart,
                plotTEnd)
        signalPlot(t, VmMiddle, ax_Vm, labelx='', labely = 'I cell $V_m$',
                nThetaTicks=5)
        signalPlot(t, VmEdge, ax_Vm, labelx='', labely = 'I cell $V_m$',
            nThetaTicks=5)
        plt.xlim(x_lim)

        # Plot E Isyn
        ax_Isyn = subplot2grid((rows, cols), (1, 3), colspan=3)
        t, IsynMiddle, IsynEdge = self._extractStateVars(mon_e, \
                ['I_clamp_GABA_A'], plotTStart, plotTEnd)
        signalPlot(t, IsynMiddle, ax_Isyn, labelx="", labely =
                'E cell $I_{\mathrm{syn}}$', nThetaTicks=5)
        signalPlot(t, IsynEdge, ax_Isyn, labelx="", labely =
                'E cell $I_{\mathrm{syn}}$', nThetaTicks=5)
        plt.xlim(x_lim)

        # Plot I Isyn
        ax_Isyn = subplot2grid((rows, cols), (2, 3), colspan=3)
        t, IsynMiddle, IsynEdge = self._extractStateVars(mon_i, \
                ['I_clamp_AMPA', 'I_clamp_NMDA'], plotTStart, plotTEnd)
        signalPlot(t, IsynMiddle, ax_Isyn, labelx='',
                labely = 'I cell $I_{syn}$', nThetaTicks=5)
        signalPlot(t, IsynEdge, ax_Isyn, labelx='',
                labely = 'I cell $I_{syn}$', nThetaTicks=5)
        plt.xlim(x_lim)


        # Check out what bump firing rate times are
        if (self.bumpTStart is None):
            self.bumpTStart = self.getOption(data, 'theta_start_t')
        if (self.bumpTEnd is None):
            self.bumpTEnd = simT

        # Plot E FRs
        ax_Isyn = subplot2grid((rows, cols), (1, 6), colspan=2)
        sp = self._getSpikeTrain(data, 'spikeMon_e', ['Ne_x', 'Ne_y'])
        Fe = sp.avgFiringRate(self.bumpTStart, self.bumpTEnd)
        Ne_x = self.getNetParam(data, 'Ne_x')
        Ne_y = self.getNetParam(data, 'Ne_y')
        bump_e = np.reshape(Fe, (Ne_y, Ne_x))
        torusFiringRate(
                rateMap  = bump_e,
                labelx   = '',
                labely   = 'Neuron #',
                titleStr = 'E firing rate')

        # Plot I FRs
        ax_Isyn = subplot2grid((rows, cols), (2, 6), colspan=2)
        sp = self._getSpikeTrain(data, 'spikeMon_i', ['Ni_x', 'Ni_y'])
        Fi = sp.avgFiringRate(self.bumpTStart, self.bumpTEnd)
        Ni_x = self.getNetParam(data, 'Ni_x')
        Ni_y = self.getNetParam(data, 'Ni_y')
        bump_i = np.reshape(Fi, (Ni_y, Ni_x))
        torusFiringRate(
                rateMap  = bump_i,
                labelx   = 'Neuron #',
                labely   = 'Neuron #',
                titleStr = 'I firing rate')

        # Plot a sample autocorrelation function
        ax_AC = subplot2grid((rows, cols), (3, 0), colspan=3)
        acVec = a['acVec']
        acdt = data['stateMonF_e'][0]['interval']
        acTimes = np.arange(acVec.shape[1]) * acdt
        signalPlot(acTimes, acVec[0], ax_AC, labelx=None,
                labely = 'I correlation')
        plt.axis('tight')
        plt.ylim([-1, 1])


        # Plot all correlations
        ax_AC = subplot2grid((rows, cols), (3, 3), colspan=3)
        acVec = a['acVec']
        acdt = data['stateMonF_e'][0]['interval']
        acTimes = np.arange(acVec.shape[1]) * acdt
        plt.hold('on')
        for nIdx in xrange(acVec.shape[0]):
            signalPlot(acTimes, acVec[nIdx], ax_AC, labelx=None, labely = '')
        plt.axis('tight')
        plt.ylim([-1, 1])


        # Plot annotations
        # TODO
        ax_ann = subplot2grid((rows, cols), (0, 0), colspan=5)
        #ax_ann.xaxis.set_visible(False)
        #ax_ann.yaxis.set_visible(False)
        ax_ann.axison = False
        txt1 = 1.00
        txt2 = 0.85
        txt3 = 0.60
        txt4 = 0.45
        txt5 = 0.20
        txt6 = 0.05
        va = 'top'; ha='left'
        ax_ann.text(0, txt1, 'Avg. freq: ',                 va=va, ha=ha)
        ax_ann.text(0, txt2, 'Avg. correlation: ',          va=va, ha=ha)
        ax_ann.text(0, txt3, 'Bump $\sigma$: ',             va=va, ha=ha)
        ax_ann.text(0, txt4, 'Bump error of fit: ',         va=va, ha=ha)
        ax_ann.text(0, txt5, 'Avg. E cell firing rate: ',   va=va, ha=ha)
        ax_ann.text(0, txt6, 'Avg. I cell firing rate: ',   va=va, ha=ha)
        
        freq      = np.mean(self._dictData(a, ['freq']))
        C         = np.mean(self._dictData(a, ['acVal']))
        bumpSigma = self._dictData(a, ['bump_e', 'sigma'])
        bumpErr   = np.sqrt(self._dictData(a, ['bump_e', 'err2']))
        FRe       = self._dictData(a, ['FR_e', 'avg'])
        FRi       = self._dictData(a, ['FR_i', 'avg'])
        ax_ann.text(0.3, txt1, "{0:.2f} Hz".format(freq),           va=va, ha=ha)
        ax_ann.text(0.3, txt2, "{0:.2f} ".format(C),                va=va, ha=ha)
        ax_ann.text(0.3, txt3, "{0:.2f} neurons".format(bumpSigma), va=va, ha=ha)
        ax_ann.text(0.3, txt4, "{0:.2f} Hz".format(bumpErr),        va=va, ha=ha)
        ax_ann.text(0.3, txt5, "{0:.2f} Hz".format(FRe),            va=va, ha=ha)
        ax_ann.text(0.3, txt6, "{0:.2f} Hz".format(FRi),            va=va, ha=ha)

        ax_ann.text(0.5, txt1, 'E total: ', va=va, ha=ha)
        ax_ann.text(0.5, txt2, 'I total: ', va=va, ha=ha)
        E_total = self.getOption(data, 'g_AMPA_total')
        I_total = self.getOption(data, 'g_GABA_total')
        ax_ann.text(0.6, txt1, "{0:.2f} nS".format(E_total), va=va, ha=ha,
                weight='bold')
        ax_ann.text(0.6, txt2, "{0:.2f} nS".format(I_total), va=va, ha=ha,
                weight='bold')

        plt.tight_layout()

        jobNum = self.getOption(data, 'job_num')
        fileName = \
            "{0}/job{1:05}_trial{2:03}_E_{3}_I_{4}_detailPlot.png".format(self.rootDir,
                    jobNum, trialNum, int(E_total), int(I_total))
        plt.savefig(fileName)



class GridPlotVisitor(DictDSVisitor):
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


