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
from matplotlib.pyplot import figure, plot, pcolormesh, subplot2grid
import os
import errno

from interface       import DictDSVisitor, sumAllVariables, \
        extractStateVariable
from plotting.signal import signalPlot
from analysis.spikes import PopulationSpikes
from plotting.bumps  import torusFiringRate
from otherpkg.log    import log_warn

__all__ = ['DetailedPlotVisitor']


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
        fig = figure(figsize=(12, 8))
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
        cols = 5
        ax_Vm = subplot2grid((rows, cols), (1, 0), colspan=2)
        t, VmMiddle, VmEdge = self._extractStateVars(mon_e, ['V_m'], plotTStart,
                plotTEnd)
        signalPlot(t, VmMiddle, ax_Vm, labelx="", labely = 'E cell $V_m$')
        signalPlot(t, VmEdge, ax_Vm, labelx="", labely = 'E cell $V_m$')
        plt.xlim(x_lim)

        # Plot I Vm
        ax_Vm = subplot2grid((rows, cols), (2, 0), colspan=2)
        t, VmMiddle, VmEdge = self._extractStateVars(mon_i, ['V_m'], plotTStart,
                plotTEnd)
        signalPlot(t, VmMiddle, ax_Vm, labelx='', labely = 'I cell $V_m$')
        signalPlot(t, VmEdge, ax_Vm, labelx='', labely = 'I cell $V_m$')
        plt.xlim(x_lim)

        # Plot E Isyn
        ax_Isyn = subplot2grid((rows, cols), (1, 2), colspan=2)
        t, IsynMiddle, IsynEdge = self._extractStateVars(mon_e, \
                ['I_clamp_GABA_A'], plotTStart, plotTEnd)
        signalPlot(t, IsynMiddle, ax_Isyn, labelx="", labely =
                'E cell $I_{\mathrm{syn}}$')
        signalPlot(t, IsynEdge, ax_Isyn, labelx="", labely =
                'E cell $I_{\mathrm{syn}}$')
        plt.xlim(x_lim)

        # Plot I Isyn
        ax_Isyn = subplot2grid((rows, cols), (2, 2), colspan=2)
        t, IsynMiddle, IsynEdge = self._extractStateVars(mon_i, \
                ['I_clamp_AMPA', 'I_clamp_NMDA'], plotTStart, plotTEnd)
        signalPlot(t, IsynMiddle, ax_Isyn, labelx='', labely = 'I cell $I_{syn}$')
        signalPlot(t, IsynEdge, ax_Isyn, labelx='', labely = 'I cell $I_{syn}$')
        plt.xlim(x_lim)


        # Check out what bump firing rate times are
        if (self.bumpTStart is None):
            self.bumpTStart = self.getOption(data, 'theta_start_t')
        if (self.bumpTEnd is None):
            self.bumpTEnd = simT

        # Plot E FRs
        ax_Isyn = subplot2grid((rows, cols), (1, 4))
        sp = self._getSpikeTrain(data, 'spikeMon_e', ['Ne_x', 'Ne_y'])
        Fe = sp.avgFiringRate(self.bumpTStart, self.bumpTEnd)
        Ne_x = self.getNetParam(data, 'Ne_x')
        Ne_y = self.getNetParam(data, 'Ne_y')
        bump_e = np.reshape(Fe, (Ne_y, Ne_x))
        torusFiringRate(
                rateMap  = bump_e,
                labelx   = 'Neuron #',
                labely   = 'Neuron #',
                titleStr = 'E firing rate')

        # Plot I FRs
        ax_Isyn = subplot2grid((rows, cols), (2, 4))
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
        ax_AC = subplot2grid((rows, cols), (3, 0), colspan=2)
        acVec = a['acVec']
        acdt = data['stateMonF_e'][0]['interval']
        acTimes = np.arange(acVec.shape[1]) * acdt
        signalPlot(acTimes, acVec[0], ax_AC, labelx=None,
                labely = 'I correlation')
        plt.axis('tight')
        plt.ylim([-1, 1])


        # Plot all correlations
        ax_AC = subplot2grid((rows, cols), (3, 2), colspan=2)
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
        ax_ann.text(0.2, txt1, "{0:.2f} Hz".format(freq),           va=va, ha=ha)
        ax_ann.text(0.2, txt2, "{0:.2f} ".format(C),                va=va, ha=ha)
        ax_ann.text(0.2, txt3, "{0:.2f} neurons".format(bumpSigma), va=va, ha=ha)
        ax_ann.text(0.2, txt4, "{0:.2f} Hz".format(bumpErr),        va=va, ha=ha)
        ax_ann.text(0.2, txt5, "{0:.2f} Hz".format(FRe),            va=va, ha=ha)
        ax_ann.text(0.2, txt6, "{0:.2f} Hz".format(FRi),            va=va, ha=ha)

        plt.tight_layout()

        jobNum = self.getOption(data, 'job_num')
        fileName = \
            "{0}/job{1:05}_trial{2:03}_detailPlot.png".format(self.rootDir,
                    jobNum, trialNum)
        plt.savefig(fileName)


