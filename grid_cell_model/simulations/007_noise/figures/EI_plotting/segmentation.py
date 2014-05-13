#
#   segmentation.py
#
#   Parameter sweeps segmentation routines.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
import numpy.ma as ma
import matplotlib.pyplot as plt

from matplotlib           import rcParams as rcp
from plotting.global_defs import globalAxesSettings
from plotting.low_level   import zeroLines
from grid_cell_model.analysis import clustering

from . import sweeps
from . import aggregate as aggr
from .base import filterData


def plotDiffHistograms(data, noise_sigmas, **kw):
    '''
    **Parameters:**

    data : a list of 2D matrices
        2D parameter sweep for all noise_sigma levels.
    noise_sigmas : sequence of floats
        Noise sigma values, must be of the same length as data
    '''
    ax              = kw.pop('ax', plt.gca())
    which           = kw.pop('which', None)
    xlabel          = kw.pop('xlabel', '')
    ylabel          = kw.pop('ylabel', 'p($\cdot$)')
    filterThreshold = kw.pop('filterThreshold', -np.infty)
    kw['linewidth'] = kw.get('linewidth', 0.01)
    kw['edgecolor'] = kw.get('edgecolor', 'white')
    kw['rwidth']    = kw.get('rwidth', 0.8)

    if (len(noise_sigmas) != len(data)):
        raise ValueError("len(noise_sigmas) != len(data)")

    #import pdb; pdb.set_trace()
    stackedData = aggr.collapseSweeps(data)
    stackedData, _ = filterData(stackedData, filterThreshold)
    diffData = np.diff(stackedData, axis=0)

    if which is None:
        whichIdx = xrange(len(noise_sigmas) - 1)
    else:
        whichIdx = [which]

    globalAxesSettings(ax)
    for idx in whichIdx:
        nan_idx = np.logical_not(np.isnan(diffData[idx, :]))
        ax.hist(diffData[idx, nan_idx], normed=True,
                alpha=1./len(noise_sigmas), **kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax



def plotDiffScatter(data, noise_sigmas, **kw):
    '''
    2D scatter plot of differences between (150-0)pA and (300-150)pA.

    **Parameters:**

    data : a list of 2D matrices
        2D parameter sweep for all noise_sigma levels.
    noise_sigmas : sequence of floats
        Noise sigma values, must be of the same length as data
    doSegmentation : bool, optional
        Whether to do segmentation or not
    '''
    ax              = kw.pop('ax', plt.gca())
    xlabel          = kw.pop('xlabel', '')
    ylabel          = kw.pop('ylabel', '')
    zerolines       = kw.pop('zerolines', True)
    doSegmentation  = kw.pop('doSegmentation', False)
    thresholds      = kw.pop('thresholds', None)
    mergeInfo       = kw.pop('mergeInfo', None)
    filterThreshold = kw.pop('filterThreshold', -np.infty)
    kw['edgecolor'] = kw.get('edgecolor', 'white')

    if (len(noise_sigmas) != len(data)):
        raise ValueError("len(noise_sigmas) != len(data)")

    stackedData = aggr.collapseSweeps(data)
    stackedData, _ = filterData(stackedData, filterThreshold)
    diffData = np.diff(stackedData, axis=0)

    if doSegmentation:
        differences = [
                diffData[0, :],
                diffData[1, :]]
        clusters = clustering.ThresholdClusters(differences, thresholds)
        if mergeInfo is not None:
            clusters = clusters.mergeClusters(mergeInfo)
        kw['c'] = clusters.assignClusters()

    globalAxesSettings(ax)
    ax.scatter(diffData[0, :], diffData[1, :], **kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.axis('scaled')

    if (zerolines):
        zeroLines(ax)

    return ax


def plotSweepSegments(spaceList, noise_sigmas, iterList, thresholds, aggrTypes,
        NTrials, **kw):
    #kw arguments
    r               = kw.pop('r', 0)
    c               = kw.pop('c', 0)
    mergeInfo       = kw.pop('mergeInfo', None)
    filterThreshold = kw.pop('filterThreshold', -np.infty)

    # Prepare data
    data = []
    for ns_idx, _ in enumerate(noise_sigmas):
        d, _, _ = aggr.aggregateType(spaceList[ns_idx], iterList, aggrTypes,
                NTrials, ignoreNaNs=True)
        data.append(d)

    # Clustering
    stackedData = aggr.collapseSweeps(data)
    stackedData, filterIdx = filterData(stackedData, filterThreshold)
    diffData = np.diff(stackedData, axis=0)
    differences = [
            diffData[0, :],
            diffData[1, :]]
    clusters = clustering.ThresholdClusters(differences, thresholds)
    if mergeInfo is not None:
        clusters = clusters.mergeClusters(mergeInfo)

    # Add another segment - filtered data
    #nClusters  = len(clusters) + 1
    clusterIdx = clusters.assignClusters()
    #clusterIdx[filterIdx.filtered] = nClusters - 1

    # Plot
    space0 = spaceList[0]
    Y, X = aggr.computeYX(space0, iterList, r=r, c=c)
    sweepSegments = np.reshape(clusterIdx, data[0].shape)
    S, ax, cax = sweeps.plot2DTrial(X, Y, sweepSegments, **kw)

    return S, ax, cax
    
