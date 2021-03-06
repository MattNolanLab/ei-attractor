'''Initialisation routines for E-I plotting.

.. currentmodule:: noisefigs.EI_plotting.base

Classes
-------

.. autosummary::
    NoiseDataSpaces
    FilteringResult

Functions
---------

.. autosummary::

    getNoiseRootDir
    getNoiseRoots
    getNoiseDataSpaces
    getDataSpace
    generateThetaSignal
    setSignalAxes
    plotStateSignal
    plotThetaSignal
    plotOneHist
    createColorbar
    filterData
    extractRateMaps
'''
from __future__ import absolute_import, print_function, division

import numpy       as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

from matplotlib.pyplot   import gca, axis, colorbar
from matplotlib.ticker   import MaxNLocator, AutoMinorLocator, LinearLocator
from matplotlib.colorbar import make_axes

from grid_cell_model.data_storage.sim_models import ei
from grid_cell_model.analysis import spikes as aspikes
from grid_cell_model.parameters.param_space import JobTrialSpace2D
from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.plotting.low_level   import xScaleBar


theta_T = 250.0 # ms
theta_dt = 1.0 # ms
theta_f = 8.0  # Hz
theta_DC = 300.0 # pA
theta_A  = 375.0 # pA
thetaMax = 1500
thetaMin = -500
thetaLim = [thetaMin, thetaMax]


def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots


def getNoiseDataSpaces(dataRoot, noise_sigmas, shape, space_cls=None):
    '''Create a list of data spaces for each of the noise sigmas.

    Parameters
    ----------
    dataRoot : str
        Path to the root of the data.
    noise_sigmas : list of numbers
        A list containing the noise levels.
    shape : A list (usually a pair) of non-negative ints
        The shape of the parameter space.
    space_cls : a class type or None
        The class which will be used to create the parameter space object for
        each noise level.
    '''
    if space_cls is None:
        space_cls = JobTrialSpace2D

    if dataRoot is None:
        return None

    roots = getNoiseRoots(dataRoot, noise_sigmas)
    ds = []
    for root in roots:
        ds.append(space_cls(shape, root))
    return ds


def getDataSpace(dataRoot, shape, space_cls=None):
    '''Create a data space.

    Parameters
    ----------
    dataRoot : str
        Path to the root of the data. If ``None``, nothing will be created.
    noise_sigmas : list of numbers
        A list containing the noise levels.
    shape : A list (usually a pair) of non-negative ints
        The shape of the parameter space.
    space_cls : a class type or None
        The class which will be used to create the parameter space object for
        each noise level.
    '''
    if space_cls is None:
        space_cls = JobTrialSpace2D

    if dataRoot is None:
        return None
    else:
        return space_cls(shape, dataRoot)


class NoiseDataSpaces(object):
    '''
    A container for the data spaces for all levels of noise.
    '''
    class Roots(object):
        def __init__(self, *args, **kw):
            la = len(args)
            if (la == 1):
                self.bump  = args[0]
                self.v     = args[0]
                self.grids = args[0]
            elif (la == 3):
                self.bump  = args[0]
                self.v     = args[1]
                self.grids = args[2]
            else:
                raise IndexError("Roots class constructor needs either 1 or"+\
                        " three arguments")
            for key, val in kw.iteritems():
                setattr(self, key, val)

    def __init__(self, roots, shape, noise_sigmas, space_cls=None):
        self.bumpGamma    = getNoiseDataSpaces(roots.bump,  noise_sigmas,
                                               shape, space_cls)
        self.v            = getNoiseDataSpaces(roots.v,     noise_sigmas,
                                               shape, space_cls)
        self.grids        = getNoiseDataSpaces(roots.grids, noise_sigmas,
                                               shape, space_cls)
        if hasattr(roots, 'constPos'):
            self.constPos = getNoiseDataSpaces(roots.constPos, noise_sigmas,
                                               shape, space_cls)

        self.conn = None
        if hasattr(roots, 'conn'):
            # Here, we only use the second dimension. The shape of the 1D
            # connection data space should match the parameter sweeps.
            if shape is not None:
                self.conn = getDataSpace(roots.conn, (1, shape[1]), space_cls)
            else:
                self.conn = None

        self.noise_sigmas = noise_sigmas


def getOption(data, optStr):
    return data['options'][optStr]


def generateThetaSignal(noise_sigma):
    t = np.arange(0, theta_T, theta_dt) * 1e-3
    phase = np.pi
    normSig = 0.5*(1. + np.cos(2*np.pi*theta_f*t + phase))
    noise = noise_sigma * np.random.randn(len(t))
    return t, theta_DC + theta_A*normSig + noise


def setSignalAxes(ax, leftSpineOn):
    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if (leftSpineOn == False):
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)

    ax.yaxis.set_major_locator(MaxNLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def plotStateSignal(ax, t, sig, leftSpineOn=True, labely="", labelyPos=-0.2,
        color='black', scaleBar=None, scaleX=0.75, scaleY=0.05, scaleText='ms',
        scaleTextSize='small'):
    setSignalAxes(ax, leftSpineOn)

    if (sig is not None):
        ax.plot(t, sig, color=color)

    #ax.set_ylabel(labely)
    ax.text(labelyPos, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)
    ax.set_xlim([t[0], t[-1]])

    if (scaleBar is not None):
        xScaleBar(scaleBar, x=scaleX, y=scaleY, ax=ax, size=scaleTextSize,
                unitsText=scaleText)


def plotThetaSignal(ax, t, theta, noise_sigma, yLabelOn, thetaLim, color='grey'):
    setSignalAxes(ax, leftSpineOn=False)
    ax.plot(t, theta, color=color)
    ax.set_ylim(thetaLim)
    txt = '$\sigma$ = ' + str(noise_sigma) + ' pA'
    ax.text(0.5, 1.1, txt,
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes,
            fontsize='large', fontweight='normal')
    ax.axhline(0.0, color='grey', linestyle=':', linewidth=0.5)
    if (yLabelOn):
        ax.text(t[-1] - 10, -50, "0 pA", ha="right", va='top', fontsize='small')


def plotOneHist(data, bins=40, normed=False, **kw):
    ax = kw.pop('ax', plt.gca())
    globalAxesSettings(ax)
    ax.hist(data, bins=bins, normed=normed, histtype='step',
            align='mid', **kw)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))


def createColorbar(ax, **kwargs):
    cbLabel     = kwargs.pop('label', '')
    #orientation = kwargs.get('orientation', 'horizontal')
    rasterized  = kwargs.pop('rasterized', None)
    labelpad    = kwargs.pop('labelpad', None)
    mappable    = kwargs.pop('mappable', None)

    cax, kwargs = make_axes(ax, **kwargs)
    globalAxesSettings(cax)
    fig = ax.figure
    cb = fig.colorbar(mappable, ax=ax, cax=cax, **kwargs)
    cb.set_label(cbLabel, labelpad=labelpad)
    cb.solids.set_rasterized(rasterized)

    return cax


class FilteringResult(object):
    def __init__(self, filteredIndexes, retainedIndexes, missing=None):
        self._filtered = filteredIndexes
        self._retained = retainedIndexes
        self._missing  = missing

    @property
    def filtered(self):
        return self._filtered

    @property
    def retained(self):
        return self._retained

    @property
    def missing(self):
        return self._missing


def filterData(stackedData, threshold):
    '''Filter some data.

    Gridness must be more than the threshold in at least one of the noise
    levels, otherwise the values will be masked.
    '''
    retainedIndexes = []
    filteredIndexes = []
    missingIndexes = []
    for dataIdx in xrange(stackedData.shape[1]):
        if np.all(stackedData.mask[:, dataIdx] == False):
            if not np.any(stackedData[:, dataIdx] > threshold):
                stackedData.mask[:, dataIdx] = True
                filteredIndexes.append(dataIdx)
            else:
                retainedIndexes.append(dataIdx)
        else:
            missingIndexes.append(dataIdx)
    return stackedData, FilteringResult(filteredIndexes, retainedIndexes,
            missingIndexes)


def extractRateMaps(ps, r, c, trialNum):
    data = ps[r][c][trialNum].data
    tStart = 0.
    win_dt = 125.      # ms
    tEnd = ei.getOption(data, 'time') - win_dt
    winLen = 250.0 # ms
    Ne_x = ei.getNetParam(data, 'Ne_x')
    Ne_y = ei.getNetParam(data, 'Ne_y')
    mon = data['spikeMon_e']
    senders, times = ei.extractSpikes(mon)
    pop = aspikes.TorusPopulationSpikes(senders, times, (Ne_x, Ne_y))
    return pop.slidingFiringRate(tStart, tEnd, win_dt, winLen)


