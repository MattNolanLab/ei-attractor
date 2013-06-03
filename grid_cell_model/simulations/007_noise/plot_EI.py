#!/usr/bin/env python
#
#   analysis_EI.py
#
#   Theta/gamma analysis using a custom "peak" method - E/I coupling parameter
#   sweep.
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
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt

import numpy.ma as ma
from matplotlib.ticker  import MaxNLocator, LinearLocator, AutoMinorLocator

from parameters           import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings, createColorbar
import logging as lg
lg.basicConfig(level=lg.WARN)
#lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

###############################################################################

def getDictData(d, keyList):
    ret = d
    for key in keyList:
        ret = ret[key]
    return ret

def aggregate2DTrial(sp, varList, iterList, trialNumList):
    if (len(iterList) != 2):
        raise ValueError("iterList must contain exactly 2 elements.")
    shape = sp.getShape()
    rows = shape[0]
    cols = shape[1]
    retVar = np.zeros(shape)
    for r in xrange(rows):
        for c in xrange(cols):
            if (len(sp[r][c]) == 0):
                retVar[r][c] = np.nan
            else:
                for trialNum in trialNumList:
                    data = sp[r][c][trialNum].data['analysis']
                    retVar[r][c] += np.mean(getDictData(data, varList))

    return retVar / len(trialNumList)

def plot2DTrial(X, Y, C, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="",
        clbarNTicks=2, xticks=True, yticks=True):

    ax = plt.gca()
    globalAxesSettings(ax)
    ax.minorticks_on()
    plt.pcolormesh(X, Y, C, vmin=vmin, vmax=vmax)
    if (clbarNTicks == None):
        createColorbar(ax, None, clBarLabel, orientation='horizontal', pad=0.2)
    else:
        createColorbar(ax, C, clBarLabel, nticks=clbarNTicks,
                orientation='horizontal', pad=0.2)
    if (xlabel != ""):
        plt.xlabel(xlabel, va='top')
        ax.xaxis.set_label_coords(0.5, -0.15)
    if (ylabel != ""):
        plt.ylabel(ylabel, ha='right')
        ax.yaxis.set_label_coords(-0.10, 0.5)
    ax.xaxis.set_ticks([0, 6000])
    ax.yaxis.set_ticks([0, 4000])
    ax.xaxis.set_minor_locator(AutoMinorLocator(6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    plt.axis('scaled')
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return C


def plotACTrial(sp, varList, iterList, trialNumList=[0], xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, iterList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = sp.getIteratedParameters(iterList)
    plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin, vmax,
            title, clbarNTicks, xticks, yticks)

def plotBumpSigmaTrial(sp, varList, iterList, thr=np.infty, trialNumList=[0],
        xlabel="", ylabel="", colorBar=True, clBarLabel="", vmin=None,
        vmax=None, title="", clbarNTicks=2, xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, iterList, trialNumList)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = sp.getIteratedParameters(iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotBumpErrTrial(sp, varList, iterList, thr=np.infty, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    C = np.sqrt(aggregate2DTrial(sp, varList, iterList, trialNumList))
    if mask is None:
        mask = False
    C = ma.MaskedArray(C, mask=np.logical_or(np.logical_or(np.isnan(C), C >
        thr), mask))
    Y, X = sp.getIteratedParameters(iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotFRTrial(sp, varList, iterList, thr=np.infty, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    FR = aggregate2DTrial(sp, varList, iterList, trialNumList)
    if mask is None:
        mask = False
    FR = ma.MaskedArray(FR, mask=np.logical_or(FR > thr, mask))
    Y, X = sp.getIteratedParameters(iterList)
    return plot2DTrial(X, Y, FR, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)
            


###############################################################################

dirs = \
    ('2013-04-24T21-43-32_EI_param_sweep_300pA_big',    (40, 40))
    #('2013-04-24T21-37-47_EI_param_sweep_150pA_big',  (40, 40))
    #('2013-04-24T15-27-30_EI_param_sweep_0pA_big',    (40, 40))
    #('2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample', (2, 2))

NTrials = 5
rootDir = "output_local/{0}".format(dirs[0])
shape   = dirs[1]

sp = JobTrialSpace2D(shape, rootDir)
iterList  = ['g_AMPA_total', 'g_GABA_total']
    
################################################################################
plt.figure(figsize=(5.1, 2.9))
N = 2
plt.subplot(1, N, 1)

acVal = plotACTrial(sp, ['acVal'], iterList,
        trialNumList=range(NTrials),
        xlabel="I (nS)",
        ylabel='E (nS)',
        clBarLabel = "Correlation",
        clbarNTicks=3,
        vmin=0,
        vmax=0.75)
###############################################################################
plt.subplot(1, N, 2)
freq = plotACTrial(sp, ['freq'], iterList,
        trialNumList=range(NTrials),
        xlabel="I (nS)",
        clBarLabel = "Frequency (Hz)",
        clbarNTicks=3,
        yticks=False,
        vmin=0,
        vmax=150)
###############################################################################
plt.tight_layout()
noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA.png'.format(int(noise_sigma)))

################################################################################
# The following is an average over trials
################################################################################
plt.figure(figsize=(5.1, 2.9))
N = 2
plt.subplot(1, N, 1)

bumpSigmaThreshold = 50
bump_sigma = plotBumpSigmaTrial(sp, ['bump_e', 'sigma'], iterList,
        thr=bumpSigmaThreshold,
        trialNumList=range(NTrials),
        xlabel="I (nS)",
        ylabel='E (nS)',
        clBarLabel = "Gaussian $\sigma$ (nrns)",
        vmin=0,
        vmax=bumpSigmaThreshold,
        clbarNTicks=4)
###############################################################################
plt.subplot(1, N, 2)
bump_err2 = plotBumpErrTrial(sp, ['bump_e', 'err2'], iterList,
        trialNumList=range(NTrials),
        mask=bump_sigma.mask,
        xlabel="I (nS)",
        clBarLabel = "Error of fit (Hz)",
        clbarNTicks=3,
        yticks=False,
        vmin=0,
        vmax=200)
###############################################################################
plt.tight_layout()
noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA_bump.png'.format(int(noise_sigma)))
###############################################################################
###############################################################################
plt.figure(figsize=(5.1, 2.9))
N = 2
plt.subplot(1, N, 1)

FR_e = plotFRTrial(sp, ['FR_e', 'avg'], iterList,
        trialNumList=range(NTrials),
        xlabel="I (nS)",
        ylabel='E (nS)',
        clBarLabel = "E cell firing rate (Hz)",
        clbarNTicks=4,
        vmin=0,
        vmax=50)
###############################################################################
plt.subplot(1, N, 2)
FR_threshold=250
FR_i = plotFRTrial(sp, ['FR_i', 'avg'], iterList,
        trialNumList=range(NTrials),
        thr=FR_threshold,
        xlabel="I (nS)",
        clBarLabel = "I cell firing rate (Hz)",
        clbarNTicks=3,
        vmin=0,
        vmax=FR_threshold,
        yticks=False)
###############################################################################
plt.tight_layout()
noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
plt.savefig(sp.rootDir +
    '/../analysis_EI_{0}pA_FR.png'.format(int(noise_sigma)))

