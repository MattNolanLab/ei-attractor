'''Raster plots for E-I parameter sweeps.

.. currentmodule:: noisefigs.EI_plotting.rasters

Functions
---------

.. autosummary::

    EIRaster
    plotEIRaster
    plotAvgFiringRate
    plot_avg_firing_rate_spikes
'''
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

import grid_cell_model.analysis.spikes as aspikes
from grid_cell_model.plotting.global_defs       import globalAxesSettings
from grid_cell_model.plotting.low_level         import xScaleBar
from grid_cell_model.data_storage.sim_models.ei import MonitoredSpikes
from grid_cell_model.plotting.signal            import signalPlot

logger = logging.getLogger(__name__)

def getSpikes(space, spaceType, r, c, trialNum, **kw):
    if (spaceType == 'velocity'):
        IvelIdx = kw.get('IvelIdx', 5) # Should correspond to Ivel==50 pA
        data = space[r][c][trialNum].data['IvelData'][IvelIdx]
    else:
        data = space[r][c][trialNum].data
    ESpikes = MonitoredSpikes(data, 'spikeMon_e', 'net_Ne')
    try:
        ISpikes = MonitoredSpikes(data, 'spikeMon_i', 'net_Ni')
    except KeyError:
        # This is just a stub: has to be fixed
        ISpikes = aspikes.PopulationSpikes(ESpikes.N, [], [])
    return ESpikes, ISpikes


def EIRaster(space, spaceType, r, c, tLimits, trialNum=0, **kw):
    '''Extract E and I spikes from the parameter space and plot them.'''
    ESpikes, ISpikes = getSpikes(space, spaceType, r, c, trialNum, **kw)
    ax = plotEIRaster(ESpikes, ISpikes, tLimits, **kw)

    return ax


def plotEIRaster(ESpikes, ISpikes, tLimits, ylabel=None, **kw):
    '''Plot a raster plot E and I spikes.'''
    # kw arguments
    ax           = kw.pop('ax', plt.gca())
    yticks       = kw.pop('yticks', True)
    yticks_style = kw.pop('yticks_style', 'separate')
    ylabelPos    = kw.pop('ylabelPos', -0.22)
    EColor       = kw.pop('ecolor', 'red')
    IColor       = kw.pop('icolor', 'blue')
    title        = kw.pop('title', True)
    ann          = kw.pop('ann', False)
    ann_EI       = kw.pop('ann_EI', False)
    sigmaTitle   = kw.pop('sigmaTitle', True)
    noise_sigma  = kw.pop('noise_sigma', None)
    scaleBar     = kw.pop('scaleBar', None)
    scaleX       = kw.pop('scaleX', 0.75)
    scaleY       = kw.pop('scaleY', 0.05)
    scaleText    = kw.pop('scaleText', 'ms')
    scaleTextYOffset = kw.pop('scaleTextYOffset', .075)
    scaleHeight  = kw.pop('scaleHeight', .02)
    scaleTextSize= kw.pop('scaleTextSize', 'small')
    reshape_senders = kw.pop('reshape_senders', True)
    kw['markersize'] = kw.get('markersize', 1.0)

    if ylabel is None:
        ylabel = 'Neuron #'

    ESpikes = ESpikes.windowed(tLimits)
    ISpikes = ISpikes.windowed(tLimits)

    ESenders, ETimes = ESpikes.rasterData()
    ISenders, ITimes = ISpikes.rasterData()

    if reshape_senders:
        # TO REMOVE: this is to transform the neuron number from row-wise to
        # column wise indexes. A better solution has to be devised in the
        # future
        logger.warn('Nx, Ny are fixed in the code. Make sure the torus size '
                    'is the same as specified here.')
        Nx, Ny = 34, 30
        N = Nx * Ny
        if (N != ESpikes.N or N != ISpikes.N):
            raise ValueError("Fix the number of neurons in plotEIRaster")
        Ex = ESenders % Nx
        Ey = ESenders // Ny
        ESenders = Ey + Ex * Ny
        Ix = ISenders % Nx
        Iy = ISenders // Ny
        ISenders = Iy + Ix * Ny

    ISenders += ESpikes.N

    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(ti.LinearLocator(2))
    ax.yaxis.set_minor_locator(ti.NullLocator())

    ax.plot(ETimes, ESenders+1, '.', color='red',  mec='none', **kw)
    ax.plot(ITimes, ISenders+1, '.', color='blue', mec='none', **kw)

    ax.set_xlim(tLimits)
    ax.set_ylim([1, ESpikes.N+ISpikes.N])
    if (yticks_style == 'separate'):
        ax.set_yticks([1, ESpikes.N, ESpikes.N+ISpikes.N])
    ax.invert_yaxis()
    ax.text(ylabelPos, 0.5, ylabel, va='center', ha='center',
            transform=ax.transAxes, rotation=90)
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    # Annotations
    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), y=1.02,
                va='bottom', ha='center')
    if (ann):
        Y, X = aggr.computeYX(space, iterList, r=r, c=c)
        gE = Y[r, c]
        gI = X[r, c]
        txt = '$g_E$ = {0} nS\n$g_I$ = {1} nS'.format(gE, gI)
        ax.text(0.99, 1.02, txt, va='bottom', ha='right', size='small',
                transform=ax.transAxes)

    if (ann_EI):
        ax.text(-0.05, 0.75, 'E', va='center', ha='center', size='small',
                transform=ax.transAxes, color='red', weight='bold')
        ax.text(-0.05, 0.25, 'I', va='center', ha='center', size='small',
                transform=ax.transAxes, color='blue', weight='bold')

    if (scaleBar is not None):
        xScaleBar(scaleBar, x=scaleX, y=scaleY, ax=ax, size=scaleTextSize,
                unitsText=scaleText, textYOffset=scaleTextYOffset,
                height=scaleHeight)


    return ax


def plotAvgFiringRate(space, spaceType, noise_sigma, popType, r, c, tLimits,
                      trialNum=0, **kw):
    '''Plot population-average firing rates of either E or I populations.'''
    ESpikes, ISpikes = getSpikes(space, spaceType, r, c, trialNum)
    if (popType == 'E'):
        spikes = ESpikes
    elif (popType == 'I'):
        spikes = ISpikes

    return plot_avg_firing_rate_spikes(spikes, tLimits, **kw)


def plot_avg_firing_rate_spikes(spikes, tLimits, **kw):
    '''Plot population-average firing rate from PopulationSpikes.'''
    # keyword arguments
    ax           = kw.pop('ax', plt.gca())
    sigmaTitle   = kw.pop('sigmaTitle', False)
    dt           = kw.pop('dt', .5)
    winLen       = kw.pop('winLen', 2.0)
    kw['xlabel'] = False
    kw['ylabel'] = kw.get('ylabel', 'r (Hz)')

    tStart = tLimits[0]
    tEnd   = tLimits[1]

    rate, times = spikes.slidingFiringRate(tStart, tEnd, dt, winLen)
    meanRate = np.mean(rate, axis=0)

    signalPlot(times, meanRate, ax, **kw)

    ax.set_ylim([0, None])
    #ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    ax.yaxis.set_major_locator(ti.LinearLocator(2))

    # Annotations
    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), y=1.02,
                va='bottom', ha='center')

    return ax


