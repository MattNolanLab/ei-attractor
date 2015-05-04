'''
Noise publication figures: connection weight figures.
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from matplotlib.gridspec import GridSpec

from grid_cell_model.parameters.param_space import JobTrialSpace2D, DataSpace
import grid_cell_model.plotting.connections as pconn
from grid_cell_model.analysis.image import Position2D
from simtools.plotting.plotters import FigurePlotter

from ..EI_plotting import aggregate as aggr

DS = DataSpace

__all__ = [
    'WeightExamplesHists',
    'WeightOutE2IPlotter',
    'WeightOutI2EPlotter',
    'WeightInE2IPlotter',
    'WeightInI2EPlotter',
    'WeightGridPlotter',
    'Burak2009ConnectionPlotter',
]


iterList  = ['g_AMPA_total', 'g_GABA_total']


def plotEToI(sp, gIdx, neuronIdx, trialNum=0, **kw):
    title = kw.pop('title', 'I cell')
    ylim  = kw.pop('ylim', None)

    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_IE']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
            title=title, **kw)
    annG = gE[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    #ann = '$g_E$ = {0} nS'.format(annG)
    #ax.text(0.95, 0.9, ann, ha='right', va='bottom', fontsize='x-small',
    #        transform=ax.transAxes)
    ax.set_xlim([0, annG])
    ax.set_xticks([0, annG])
    ax.xaxis.set_ticklabels([0, '$g_E$'])
    ax.set_ylim(ylim)


def plotIToE(sp, gIdx, neuronIdx, trialNum=0, **kw):
    title = kw.pop('title', 'E cell')

    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
            title=title, **kw)
    annG = gI[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_I$ = {0} nS'.format(annG)
    ax.text(0.95, 0.9, ann, ha='right', va='bottom', fontsize='x-small',
            transform=ax.transAxes)
    ax.set_xlim([0, annG])
    ax.set_xticks([0, annG])


def plotIToEBrokenAxis(sp, gIdx, neuronIdx, trialNum=0, axBoundaries=None,
                       axesProportions=(0.5, 0.5), bottomLimits=None,
                       topLimits=None, **kw):
    if axBoundaries is None:
        axBoundaries = [0, 0, 1, 1]
    left, bottom, right, top = axBoundaries
    title = kw.pop('title', 'E cell')
    fig   = kw.pop('fig', plt.gcf())
    h = top - bottom
    w = right - left
    hBottom = h*axesProportions[0]
    hTop = h*axesProportions[1]

    axBottom = fig.add_axes(Bbox.from_extents(left, bottom, right, bottom +
                                              hBottom))
    axTop = fig.add_axes(Bbox.from_extents(left, top - hTop, right, top),
                         sharex=axBottom)

    _, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]

    pconn.plotConnHistogram(conns, title=title, ax=axBottom, **kw)
    kw['ylabel'] = ''
    pconn.plotConnHistogram(conns, title=title, ax=axTop, **kw)
    annG = gI[0, gIdx]
    if annG - int(annG) == 0:
        annG = int(annG)
    #ann = '$g_I$ = {0} nS'.format(annG)
    #fig.text(left+0.95*w, bottom+0.9*h, ann, ha='right', va='bottom',
    #        fontsize='x-small')

    axBottom.set_xlim([0, annG])
    axBottom.set_xticks([0, annG])
    axBottom.xaxis.set_ticklabels([0, '$g_I$'])
    axBottom.set_ylim(bottomLimits)
    axBottom.set_yticks(bottomLimits)
    axBottom.yaxis.set_minor_locator(ti.NullLocator())
    axTop.set_ylim(topLimits)
    axTop.set_yticks([topLimits[1]])
    axTop.xaxis.set_visible(False)
    axTop.spines['bottom'].set_visible(False)

    divLen = 0.07
    d = .015
    kwargs = dict(transform=fig.transFigure, color='k', clip_on=False)
    axBottom.plot((left-divLen*w, left+divLen*w), (bottom+hBottom + d,
                                                   bottom+hBottom - d),
                  **kwargs)
    axTop.plot((left-divLen*w, left+divLen*w), (top-hTop + d, top-hTop - d),
               **kwargs)

    return axBottom, axTop


class WeightExamplesHists(FigurePlotter):
    left   = 0.35
    bottom = 0.32
    right  = 0.95
    top    = 0.85

    neuronIdx = 0

    def __init__(self, *args, **kwargs):
        super(WeightExamplesHists, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']

        exampleFigSize = (1.6, 1.4)
        exLeft   = 0.4
        exBottom = 0.32
        exRight  = 0.95
        exTop    = 0.85
        exampleRC = ( (5, 15), (15, 5) )
        for exIdx, example in enumerate(exampleRC):
            kw = dict()
            if exIdx == 1:
                kw['xlabel'] = ''

            fig = self._get_final_fig(exampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(exLeft, exBottom, exRight, exTop))
            plotEToI(ps.conn, example[0], self.neuronIdx, ylabel='', title='',
                    rwidth=0.8,
                    linewidth=0,
                    **kw)
            ax.yaxis.set_minor_locator(ti.NullLocator())
            ax.set_xlabel(ax.xaxis.get_label_text(), labelpad=-5)
            fname = output_dir + "/figure_connections_examples_E2I{0}.pdf"
            plt.savefig(fname.format(exIdx), dpi=300, transparent=True)
            plt.close()


            fig = self._get_final_fig(exampleFigSize)
            axBoundaries = (exLeft, exBottom, exRight, exTop)
            axBottom, axTop = plotIToEBrokenAxis(ps.conn, example[1],
                                                 self.neuronIdx,
                                                 ylabel='', title='',
                                                 axBoundaries=axBoundaries,
                                                 axesProportions=(0.75, 0.2),
                                                 bottomLimits=(0, 60),
                                                 topLimits=(800, 900),
                                                 rwidth=0.8,
                                                 linewidth=0,
                                                 **kw)
            axBottom.set_xlabel(axBottom.xaxis.get_label_text(), labelpad=-5)
            fig.text(exLeft - 0.27, 0.5*(self.bottom+self.top), 'Count',
                    rotation=90, ha='center', va='center')
            fname = output_dir + "/figure_connections_examples_I2E{0}.pdf"
            plt.savefig(fname.format(exIdx), dpi=300, transparent=True)
            plt.close()


class WeightPlotter(FigurePlotter):
    '''Color plots of incoming/outgoing synaptic weights of a single neuron.'''
    def __init__(self, *args, **kwargs):
        super(WeightPlotter, self).__init__(*args, **kwargs)

    def _get_plotting_kwargs(self):
        return {
            'xlabel': self.myc.get('xlabel', 'Neuron #'),
            'ylabel': self.myc.get('ylabel', 'Neuron #'),
            'use_title': self.myc.get('use_title', True),
        }

    def plotOutgoing(self, gIdx, type, neuronIdx, trialNum=0, **kw):
        '''Plot outgoing weights from a single neuron to all other neurons.'''
        Nx = None
        Ny = None
        use_title = kw.pop('use_title', True)

        data = self.env.ps.conn[0][gIdx][trialNum].data
        if type == 'E':
            var = 'g_IE'
            Nx = DS.getNetParam(data, 'Ni_x')
            Ny = DS.getNetParam(data, 'Ni_y')
            kw['title'] = 'E cell $\\rightarrow$ I cells'

        elif type == 'I':
            var = 'g_EI'
            Nx = DS.getNetParam(data, 'Ne_x')
            Ny = DS.getNetParam(data, 'Ne_y')
            kw['title'] = 'I cell $\\rightarrow$ E cells'

        if not use_title:
            kw['title'] = ''

        conns = np.reshape(data[var][:, neuronIdx], (Ny, Nx))
        pconn.plot2DWeightMatrix(conns, **kw)

    def plotIncoming(self, gIdx, type, neuronIdx, trialNum=0, **kw):
        '''Plot incoming weights of a single neuron from all other neurons.'''
        Nx = None
        Ny = None
        use_title = kw.pop('use_title', True)

        data = self.env.ps.conn[0][gIdx][trialNum].data
        if type == 'I':
            var = 'g_IE'
            Nx = DS.getNetParam(data, 'Ne_x')
            Ny = DS.getNetParam(data, 'Ne_y')
            kw['title'] = 'E cells $\\rightarrow$ I cell'

        elif type == 'E':
            var = 'g_EI'
            Nx = DS.getNetParam(data, 'Ni_x')
            Ny = DS.getNetParam(data, 'Ni_y')
            kw['title'] = 'I cells $\\rightarrow$ E cell'

        if not use_title:
            kw['title'] = ''

        conns = np.reshape(data[var][neuronIdx, :], (Ny, Nx))
        pconn.plot2DWeightMatrix(conns, **kw)


class WeightOutE2IPlotter(WeightPlotter):
    '''Outgoing E-->I connections.'''
    def __init__(self, *args, **kwargs):
        super(WeightOutE2IPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig_size = self.myc['fig_size']
        g_idx = self.myc['g_idx']
        neuron_idx = self.myc['neuron_idx']

        # E-->I
        fig = plt.figure(figsize=fig_size)
        self.plotOutgoing(g_idx, "E", neuron_idx,
                          **self._get_plotting_kwargs())
        fig.tight_layout()
        fname = self.get_fname("/connections_pcolor_out_E2I.pdf")
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()


class WeightOutI2EPlotter(WeightPlotter):
    '''Outgoing I-->E connections.'''
    def __init__(self, *args, **kwargs):
        super(WeightOutI2EPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig_size = self.myc['fig_size']
        g_idx = self.myc['g_idx']
        neuron_idx = self.myc['neuron_idx']

        # I-->E
        fig = plt.figure(figsize=fig_size)
        self.plotOutgoing(g_idx, "I", neuron_idx,
                          **self._get_plotting_kwargs())
        fig.tight_layout()
        fname = self.get_fname("/connections_pcolor_out_I2E.pdf")
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()


class WeightInE2IPlotter(WeightPlotter):
    '''Incoming E-->I connections.'''
    def __init__(self, *args, **kwargs):
        super(WeightInE2IPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig_size = self.myc['fig_size']
        g_idx = self.myc['g_idx']
        neuron_idx = self.myc['neuron_idx']

        # Out of curiosity: plot the weights to one neuron (incoming)
        # E-->I
        fig = plt.figure('g_in_E2I', figsize=fig_size)
        self.plotIncoming(g_idx, "I", neuron_idx,
                          **self._get_plotting_kwargs())
        fig.tight_layout()
        fname = self.get_fname("/connections_pcolor_in_E2I.pdf")
        plt.savefig(fname, dpi=300, transparent=True)


class WeightInI2EPlotter(WeightPlotter):
    '''Incoming I-->E connections.'''
    def __init__(self, *args, **kwargs):
        super(WeightInI2EPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig_size = self.myc['fig_size']
        g_idx = self.myc['g_idx']
        neuron_idx = self.myc['neuron_idx']

        # I-->E
        fig = plt.figure('g_in_I2E', figsize=fig_size)
        self.plotIncoming(g_idx, "E", neuron_idx,
                          **self._get_plotting_kwargs())
        fig.tight_layout()
        fname = self.get_fname("/connections_pcolor_in_I2E.pdf")
        plt.savefig(fname, dpi=300, transparent=True)


class WeightGridPlotter(WeightPlotter):
    '''Color plots of weights of selected neurons in a grid.'''
    def __init__(self, *args, **kwargs):
        super(WeightGridPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        g_idx = self.myc['g_idx']
        neuron_idx = self.myc['neuron_idx']
        l, b, r, t = self.myc['bbox_rect']

        fig = self._get_final_fig(self.myc['fig_size'])
        gs = GridSpec(2, 2)
        gs.update(left=l, right=r, bottom=b, top=t, hspace=0)

        # E-->I outgoing
        ax = fig.add_subplot(gs[0, 0])
        self.plotOutgoing(g_idx, "E", neuron_idx, ax=ax, xlabel='', ylabel='',
                          use_title=False)
        ax.set_xticks([])

        # I-->E input
        ax = fig.add_subplot(gs[0, 1])
        self.plotIncoming(g_idx, "E", neuron_idx, ax=ax, ylabel='', xlabel='',
                          use_title=False)
        ax.set_xticks([])
        ax.set_yticks([])

        # I-->E outgoing
        ax = fig.add_subplot(gs[1, 0])
        self.plotOutgoing(g_idx, "I", neuron_idx, ax=ax, use_title=False,
                          xlabel='', ylabel='')

        # E-->I input
        ax = fig.add_subplot(gs[1, 1])
        self.plotIncoming(g_idx, "I", neuron_idx, ax=ax, xlabel='', ylabel='',
                          use_title=False)
        ax.set_yticks([])

        fname = self.get_fname("/connections_pcolor_grid.pdf")
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()

        # Add an extra colorbar
        fig = self._get_final_fig(self.myc['cbar_fig_size'])
        ax_cbar = fig.add_axes([0.05, 0.80, 0.8, 0.15])
        cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=mpl.cm.jet,
                                         norm=mpl.colors.Normalize(vmin=0,
                                                                   vmax=1),
                                         ticks=[0, 1],
                                         orientation='horizontal')
        ax_cbar.xaxis.set_ticklabels(['0', '$g_{E/I}$'])
        fname_cbar = self.get_fname("/connections_pcolor_grid_colorbar.pdf")
        plt.savefig(fname_cbar, dpi=300, transparent=True)
        plt.close()


class Burak2009ConnectionPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(Burak2009ConnectionPlotter, self).__init__(*args, **kwargs)

    def compute_weights(self, X1, X2, a, gamma, beta, l, pref_theta):
        '''Compute the outgoing weights between neurons at positions specified
        by ``X1`` and ``X2``. See Burak and Fiete (2009).
        '''
        X_sq = ((X1.x - X2.x - l * pref_theta.x)**2 +
                (X1.y - X2.y - l * pref_theta.y)**2)
        return a * np.exp(-gamma * X_sq) - np.exp(-beta * X_sq)

    def plot(self, *args, **kwargs):
        output_dir = self.config['output_dir']

        lambda_net = 20.
        a = 1.
        beta = 3. / lambda_net**2
        gamma = 1.05 * beta
        l = 5.
        X1 = Position2D(0., 0.)
        n_range = 30
        X2_x, X2_y = np.meshgrid(np.arange(-n_range, n_range),
                                 np.arange(-n_range, n_range))
        X2 = Position2D(X2_x, X2_y)

        fig = self._get_final_fig(self.myc['fig_size'])
        # Shift up
        ax_up = fig.add_subplot(2, 2, 1)
        pref_theta = Position2D(0, -1)
        w = self.compute_weights(X1, X2, a, gamma, beta, l, pref_theta)
        ax_up.pcolor(X2.x, X2.y, w, rasterized=True)
        ax_up.set_xticks([])
        ax_up.set_yticks([])


        # Shift down
        ax_down = fig.add_subplot(2, 2, 2)
        pref_theta = Position2D(0, 1)
        w = self.compute_weights(X1, X2, a, gamma, beta, l, pref_theta)
        ax_down.pcolor(X2.x, X2.y, w, rasterized=True)
        ax_down.set_xticks([])
        ax_down.set_yticks([])

        # Shift left
        ax_left = fig.add_subplot(2, 2, 3)
        pref_theta = Position2D(1, 0)
        w = self.compute_weights(X1, X2, a, gamma, beta, l, pref_theta)
        ax_left.pcolor(X2.x, X2.y, w, rasterized=True)
        ax_left.set_xticks([])
        ax_left.set_yticks([])

        # Shift right
        ax_right = fig.add_subplot(2, 2, 4)
        pref_theta = Position2D(-1, 0)
        w = self.compute_weights(X1, X2, a, gamma, beta, l, pref_theta)
        ax_right.pcolor(X2.x, X2.y, w, rasterized=True)
        ax_right.set_xticks([])
        ax_right.set_yticks([])

        fig.tight_layout()
        fname = self.config['output_dir'] + "/intro_burak2009_conn_weights.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close()
