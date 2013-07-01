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

from parameters           import JobTrialSpace2D, DataSpace
from plotting.global_defs import globalAxesSettings, createColorbar
import logging as lg
lg.basicConfig(level=lg.WARN)
#lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

###############################################################################

def aggregate2DTrial(sp, varList, trialNumList):
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean)
    return np.mean(retVar, 2)



def computeYX(sp, iterList):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ne')
    Ni = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ni')
    return E/Ne, I/Ni

def plotACTrial(sp, varList, iterList, trialNumList=[0], xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = computeYX(sp, iterList)
    plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin, vmax,
            title, clbarNTicks, xticks, yticks)

def plotBumpSigmaTrial(sp, varList, iterList, thr=np.infty, trialNumList=[0],
        xlabel="", ylabel="", colorBar=True, clBarLabel="", vmin=None,
        vmax=None, title="", clbarNTicks=2, xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = computeYX(sp, iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotBumpErrTrial(sp, varList, iterList, thr=np.infty, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    C = np.sqrt(aggregate2DTrial(sp, varList, trialNumList))
    if mask is None:
        mask = False
    C = ma.MaskedArray(C, mask=np.logical_or(np.logical_or(np.isnan(C), C >
        thr), mask))
    Y, X = computeYX(sp, iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotFRTrial(sp, varList, iterList, thr=np.infty, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    FR = aggregate2DTrial(sp, varList, trialNumList)
    if mask is None:
        mask = False
    FR = ma.MaskedArray(FR, mask=np.logical_or(FR > thr, mask))
    Y, X = computeYX(sp, iterList)
    return plot2DTrial(X, Y, FR, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)
            

###############################################################################

dirs = \
    ('2013-04-24T15-27-30_EI_param_sweep_0pA_big',    (40, 40))
    #('2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample', (2, 2))
    #('2013-04-24T21-43-32_EI_param_sweep_300pA_big',    (40, 40))
    #('2013-04-24T21-37-47_EI_param_sweep_150pA_big',  (40, 40))

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

