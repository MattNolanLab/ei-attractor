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
from interface       import DictDSVisitor
from plotting        import signalPlot
from analysis.spikes import PopulationSpikes
from plotting.bumps  import torusFiringRate

import matplotlib.pyplot as plt
from plt import figure, plot, pcolormesh, subplot2grid

__all__ = []


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
        '''
        self.rootDir    = rootDir
        self.plotT      = plotT
        self.bumpTStart = bumpTStart
        self.bumpTEnd   = bumpTEnd

    def _getSpikeTrain(self, data, monName, dimList):
        senders, times, N = DictDSVisitor._getSpikeTrain(self, data, monName,
                dimList)
        return PopulationSpikes(N, senders, times)
    
    def visitDictDataSet(self, ds):
        data = ds.data
        fig = figure()
        mon_e = data['stateMon_e']
        mon_i = data['stateMon_i']
        nIdxMiddle = 0
        nIdxEdge = 1
        
        # Plot E Vm
        ax_Vm = subplot2grid((3, 3), (1, 0))
        t = extractStateVariable(mon_e, nIdxMiddle, 'times')
        Vm = extractStateVariable(mon_e, nIdxMiddle, 'V_m')
        signalPlot(t, Vm, ax_Vm, labelx="", labely = '$V_m$')
        Vm = extractStateVariable(mon_e, nIdxEdge, 'V_m')
        signalPlot(t, Vm, ax_Vm, labelx="", labely = '$V_m$')

        # Plot I Vm
        ax_Vm = subplot2grid((3, 3), (2, 0))
        t = extractStateVariable(mon_i, nIdxMiddle, 'times')
        Vm = extractStateVariable(mon_i, nIdxMiddle, 'V_m')
        signalPlot(t, Vm, ax_Vm, labelx=None, labely = '$V_m$')
        Vm = extractStateVariable(mon_i, nIdxEdge, 'V_m')
        signalPlot(t, Vm, ax_Vm, labelx=None, labely = '$V_m$')

        # Plot E Isyn
        ax_Isyn = subplot2grid((3, 3), (1, 1))
        t = extractStateVariable(mon_e, nIdxMiddle, 'times')
        Isyn = extractStateVariable(mon_e, nIdxMiddle, 'I_clamp_GABA_A')
        signalPlot(t, Isyn, ax_Isyn, labelx="", labely = '$I_{syn}$')
        Isyn = extractStateVariable(mon_e, nIdxEdge, 'I_clamp_GABA_A')
        signalPlot(t, Isyn, ax_Isyn, labelx="", labely = '$I_{syn}$')

        # Plot I Isyn
        ax_Isyn = subplot2grid((3, 3), (2, 1))
        t = extractStateVariable(mon_i, nIdxMiddle, 'times')
        Isyn = sumAllVariables(mon_i, nIdxMiddle, ['I_clamp_AMPA', \
                'I_clamp_NMDA'])
        signalPlot(t, Isyn, ax_Isyn, labelx=None, labely = '$I_{syn}$')
        Isyn = sumAllVariables(mon_i, nIdxEdge, ['I_clamp_AMPA', \
                'I_clamp_NMDA'])
        signalPlot(t, Isyn, ax_Isyn, labelx=None, labely = '$I_{syn}$')

        # Plot E FRs
        ax_Isyn = subplot2grid((3, 3), (1, 2))
        sp = self._getSpikeTrain(data, 'spikeMon_e', ['Ne_x', 'Ne_y'])
        Fe = sp.avgFiringRate(self.tStart, self.tEnd)
        Ne_x = self.getOption(data, 'Ne_x')
        Ne_y = self.getOption(data, 'Ne_y')
        bump_e = np.reshape(Fe, (Ne_y, Ne_x))
        torusFiringRate(
                rateMap  = bump_e,
                labelx   = 'E neuron #',
                labely   = 'I neuron #',
                titleStr = '')

        # Plot I FRs
        ax_Isyn = subplot2grid((3, 3), (2, 2))
        sp = self._getSpikeTrain(data, 'spikeMon_i', ['Ni_x', 'Ni_y'])
        Fi = sp.avgFiringRate(self.tStart, self.tEnd)
        Ni_x = self.getOption(data, 'Ni_x')
        Ni_y = self.getOption(data, 'Ni_y')
        bump_i = np.reshape(Fi, (Ni_y, Ni_x))
        torusFiringRate(
                rateMap  = bump_i,
                labelx   = 'E neuron #',
                labely   = 'I neuron #',
                titleStr = '')

        # Plot annotations
        # TODO



