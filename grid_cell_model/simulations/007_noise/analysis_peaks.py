#!/usr/bin/env python
#
#   analysis_peaks.py
#
#   Theta/gamma analysis using a custom "peak" method.
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

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from matplotlib.ticker  import MaxNLocator, LinearLocator

from plotting.global_defs import globalAxesSettings
from parameters           import JobTrialSpace2D
from analysis.visitors    import AutoCorrelationVisitor
import logging as lg
lg.basicConfig(level=lg.DEBUG)


# Other
plt.rcParams['font.size'] = 14


###############################################################################

def peakAmp(stateMon, startT, thetaFreq, fieldStr, extFunc):
    N = len(stateMon)
    thetaT = 1. / (thetaFreq * 1e-3)

    peaks = []
    for n_id in range(N):
        sig = stateMon[n_id]['events'][fieldStr]
        dt = stateMon[n_id]['interval'] # sec
        sig_split = splitSigToThetaCycles(sig[startT/dt:], thetaT, dt)
        peaks.append(globalExtremum(sig_split, extFunc))

    return np.array(peaks)


def plotNoiseSigma(noise_sigma, res_means, res_stds, newFigure=True, ylabel="",
        xlabel=None, title=""):
    if (newFigure):
        plt.figure()
    ax = plt.gca()

    X = np.arange(len(noise_sigma))
    plt.errorbar(X, res_means, res_stds, fmt='o-')
    ax.xaxis.set_ticks(X)
    ax.xaxis.set_ticklabels(np.array(noise_sigma, dtype=int))
    globalAxesSettings(ax)
    plt.margins(0.05)
    if (xlabel is None):
        plt.xlabel('$\sigma_{noise}$ (pA)')
    elif (xlabel != ""):
        plt.xlabel(xlabel)

    if (ylabel != ""):
        ax.set_ylabel(ylabel, multialignment='center')

    if (title != ""):
        plt.title(title, x=-0.3, y=1.1, ha='left', va='bottom', weight='bold',
                size='x-large')


def plotAC(AC, times, title="", TUnits='ms', ylabel="Correlation"):
    AC = np.array(AC)

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.hold('on')

    ACMean = np.mean(AC, axis=0)
    ACstd  = np.std(AC, axis=0)
    plt.fill_between(times, ACMean - ACstd, ACMean + ACstd, alpha=0.1,
            linewidth=0)
    plt.plot(times, ACMean)
    plt.xlabel('Time ({0:s})'.format(TUnits))
    if (ylabel != ""):
        plt.ylabel(ylabel)
        ytickLabel = True
    else:
        ytickLabel = False

    m = 0.01
    ax.yaxis.set_ticks([-1, 0, 1])
    if (not ytickLabel):
        ax.yaxis.set_ticklabels([])
    #plt.margins(m)
    ymin, ymax = absoluteLimits(-1, 1, m)
    plt.ylim([ymin, ymax])
    ax.xaxis.set_major_locator(LinearLocator(3))
    ax.yaxis.grid(True)

    if (title != ""):
        plt.title(title)



def plotACs(sp):
    acFigSize = (6.2, 2.25)
    plt.figure(figsize=acFigSize)
    noise_sigma = sp.getIteratedParameters(['noise_sigma'])[0].flat
    N = len(sp[0])
    
    for idx in range(N):
        acVec = getAnalysisField(sp, idx, 'acVec')
        times = np.arange(acVec.shape[1])
        if (idx == 0):
            ylabel = "Correlation"
        else:
            ylabel = ""
        plt.subplot(1, N, idx+1)
        title = '{0} pA'.format(noise_sigma[idx])
        plotAC(acVec, times, ylabel=ylabel, title=title)

    plt.tight_layout()
    plt.savefig(sp.rootDir + "/peak_analysis_AC.pdf")



def getAnalysisField(sp, idx, fieldStr, trialIdx=0):
    return sp[0][idx][trialIdx].data['analysis'][fieldStr]


def plotSigmas(sp):
    sigmaFigSize = (7,2.5)
    plt.figure(figsize=sigmaFigSize)

    noise_sigma = sp.getIteratedParameters(['noise_sigma'])[0].flat
    N = len(sp[0])
    freq_mean  = np.ndarray((N, ))
    freq_std   = np.ndarray((N, ))
    acval_mean = np.ndarray((N, ))
    acval_std  = np.ndarray((N, ))
    for idx in range(N):
        freq_mean[idx]  = np.mean(getAnalysisField(sp, idx, 'freq'))
        freq_std[idx]   = np.std(getAnalysisField(sp, idx, 'freq'))
        acval_mean[idx] = np.mean(getAnalysisField(sp, idx, 'acVal'))
        acval_std[idx]  = np.std(getAnalysisField(sp, idx, 'acVal'))
        
    plt.subplot(1, 2, 1)
    plotNoiseSigma(noise_sigma, freq_mean, freq_std,
            ylabel="Frequency (Hz)",
            title="",
            newFigure=False)
    loc = MaxNLocator(nbins=4, steps=[10])
    plt.gca().yaxis.set_major_locator(loc)
    
    plt.subplot(1, 2, 2)
    plotNoiseSigma(noise_sigma, acval_mean, acval_std, 
            ylabel="Mean\nauto-correlation",
            title="",
            newFigure=False)
    plt.tight_layout(w_pad=1., h_pad=2.5)
    
    plt.savefig(sp.rootDir + "/peak_analysis.pdf")



def computeMinMaxMargin(sig_mean, sig_std=0, margin=0):
    sm = np.array(sig_mean)
    ss  = np.array(sig_std)
    ymin, ymax = np.min(sm - ss), np.max(sm + ss)
    tot = ymax - ymin
    return ymin - tot*margin, ymax + tot*margin


def absoluteLimits(min, max, margin):
    tot = max-min
    return min-margin*tot, max+margin*tot







###############################################################################
rootDir = 'output_local/correlation_1d'
shape = (1, 3)
sp = JobTrialSpace2D(shape, rootDir)

monName   = 'stateMonF_e'
stateList = ['I_clamp_GABA_A']
iterList  = ['noise_sigma']
forceUpdate = False
visitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=forceUpdate)
sp.visit(visitor)

################################################################################

plotSigmas(sp)

plotACs(sp)

################################################################################
#plt.show()
