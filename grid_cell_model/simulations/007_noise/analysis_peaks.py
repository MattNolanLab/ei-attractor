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

from data_storage           import DataStorage
from analysis.signal        import splitSigToThetaCycles, globalExtremum, corr,\
        localExtrema, butterBandPass, autoCorrelation, downSample
from plotting.global_defs   import globalAxesSettings

import matplotlib
print matplotlib.__version__, matplotlib.__path__


# Other
plt.rcParams['font.size'] = 12


###############################################################################


if (len(sys.argv) == 1):
    raise Exception("Must specify job numbers as command line parameters!")

jobNums = []
for jn in sys.argv[1:]:
    jobNums.append(int(jn))

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

def writeArrayTxt(fileName, arr):
    f = open(fileName, 'w')
    for val in arr:
        f.write('{0:f}\n'.format(val))
    f.close()


def findFreq(ac, dt, ext_idx, ext_t):
    max_idx = np.nonzero(ext_t > 0)[0]
    if (len(max_idx) == 0):
        raise ValueError("Autocorrelation must contain at leas one local maximum")
    
    # First local maximum ac[0] excluded
    max1_idx = ext_idx[max_idx[0]]
    max1_t   = max1_idx * dt
    max1     = ac[max1_idx]

    return (1./max1_t, max1)


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
        plt.title(title, x=-0.5, y=1.1, ha='left', va='bottom', weight='bold',
                size='x-large')


def plotAC(AC, times, NExamples=None, title="", TUnits='ms',
        ylabel="Correlation"):
    if (NExamples is None):
        N = len(AC)
    elif (NExamples >= 0):
        N = NExamples
    else:
        raise ValueError("NExamples must be None or non-negative.")
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



def plotACs(AC_all, times, noise_sigma):
    acFigSize = (6.2, 2.25)
    plt.figure(figsize=acFigSize)
    N = len(AC_all)
    
    for idx in range(N):
        if (idx == 0):
            ylabel = "Correlation"
        else:
            ylabel = ""
        plt.subplot(1, N, idx+1)
        title = '{0} pA'.format(noise_sigma[idx])
        plotAC(AC_all[idx], times, ylabel=ylabel, title=title)

    plt.tight_layout()
    plt.savefig(output_dir + "/peak_analysis_AC.pdf")


def plotSigmas(noise_sigma, peaks_mean, peaks_std, freq_mean, freq_std,
        acval_mean, acval_std):
    sigmaFigSize = (7,5)
    plt.figure(figsize=sigmaFigSize)
    
    plt.subplot2grid((2, 2), (0, 0))
    plotNoiseSigma(noise_sigma, peaks_mean, peaks_std,
            xlabel = "",
            ylabel="E cell max. $I_{syn}$ / $\\theta$ \n (pA)",
            title="A",
            newFigure=False)
    plt.gca().set_xticklabels([])
    loc = MaxNLocator(nbins=4)
    plt.gca().yaxis.set_major_locator(loc)
    #ymin,ymax = computeMinMaxMargin(peaks_mean, peaks_std, margin=0.1)
    ymin,ymax = absoluteLimits(0, 11000, margin=0.05)
    plt.ylim([ymin, ymax])
    
    
    plt.subplot2grid((2, 2), (1, 0))
    plotNoiseSigma(noise_sigma, freq_mean, freq_std,
            ylabel="Frequency (Hz)",
            title="B",
            newFigure=False)
    loc = MaxNLocator(nbins=4, steps=[10])
    plt.gca().yaxis.set_major_locator(loc)
    
    plt.subplot2grid((2, 2), (1, 1))
    plotNoiseSigma(noise_sigma, acval_mean, acval_std, 
            ylabel="Mean\nauto-correlation",
            title="C",
            newFigure=False)
    plt.tight_layout(w_pad=0.0, h_pad=2.5)
    
    plt.savefig(output_dir + "/peak_analysis.pdf")



def computeMinMaxMargin(sig_mean, sig_std=0, margin=0):
    sm = np.array(sig_mean)
    ss  = np.array(sig_std)
    ymin, ymax = np.min(sm - ss), np.max(sm + ss)
    tot = ymax - ymin
    return ymin - tot*margin, ymax + tot*margin


def absoluteLimits(min, max, margin):
    tot = max-min
    return min-margin*tot, max+margin*tot



def extractStateVariable(mon, nIdx, varStr):
    '''Extract state variable from a monitor.
    
    Parameters
    ----------
    mon : dict
        A list of (NEST) monitors, each monitoring one neuron.
    nIdx : int
        Neuron index
    varStr : str
        Name of the variable
    output
        A tuple (data, dt), for the signal
    '''
    n = mon[nIdx]
    return n['events'][varStr], n['interval']


def sumAllVariables(mon, nIdx, varList):
    '''
    Extract all variables from the list of monitors and sum them
    '''
    sigSum = None
    dtCheck = None
    for idx in range(len(varList)):
        sig, dt = extractStateVariable(mon, nIdx, varList[idx])
        if (idx == 0):
            sigSum = sig
            dtCheck = dt
        else:
            assert(dtCheck == dt)
            sigSum += sig

    return sigSum, dt


def extractACStat(mon, stateList, maxLag, dtMult=1e-3, norm=True,
        bandStart=20, bandEnd=200):
    '''
    Extrac autocorrelation statistics from a monitor

    Parameters
    ----------
    mon : list of dicts
        A list of (NEST) state monitors' status dictionary
    stateList : list of strings
        A list of strings naming the state variables to extract (and sum)
    maxLag : float
        Maximal lag to extract (in seconds, NOT timesteps)
    dtMult : float, optional
        dt Multiplier to transform dt into seconds
    norm : bool, optional
        Whether the autocorrelation function should be normalized
    bandStart : float, optional
        Bandpass start frequency
    bandEnd   : float, optional
        Bandpass end frequency
    '''
    freq   = [] # Frequency of input signal
    acval  = [] # Auto-correlation at the corresponding frequency
    acVec  = []
    #for n_id in range(len(mon)):
    for n_id in range(5):
        print "n_id: ", n_id
        sig, dt = sumAllVariables(mon, n_id, stateList)
        sig = butterBandPass(sig, dt*dtMult, bandStart, bandEnd)
        ac = autoCorrelation(sig - np.mean(sig), max_lag=maxLag/dt, norm=norm)
        ext_idx, ext_t = localExtrema(ac)
        acVec.append(ac)

        #plt.figure()
        #plt.plot(times[0:slice], ac[0:slice])
        #plt.hold('on')
        #plt.plot(times[ext_idx], ac[ext_idx], '.')
        #plt.ylim([-1, 1])
        #plt.savefig(output_fname + "_peaks_ac_extrema_%d.pdf" % n_id)
        #plt.close()

        f, a = findFreq(ac, dt*1e-3, ext_idx, ext_t)
        freq.append(f)
        acval.append(a)

    return freq, acval, acVec


###############################################################################


output_dir = 'output_local/tmp_bump_fitting'
fileNamePrefix = ''

AC_all      = []
peaks_mean  = []
peaks_std   = []
noise_sigma = []
freq_mean   = []
freq_std    = []
acval_mean  = []
acval_std   = []

for jobNum in jobNums:
    print "jobNum: ", jobNum
    in_fname = "{0}/{1}job{2:05}_output".format(output_dir, fileNamePrefix, jobNum)
    output_fname = in_fname
    ds = DataStorage.open(in_fname + ".h5", 'r')
    d  = ds['trials'][0]
    options     = d['options']
    stateMonF_e = d['stateMonF_e']
    
    
    peaks = peakAmp(
            stateMon = d['stateMonF_e'],
            startT   = options['theta_start_t'],
            thetaFreq = options['theta_freq'],
            fieldStr = 'I_clamp_GABA_A',
            extFunc = np.max).flat
    peaks_mean.append(np.mean(peaks))
    peaks_std.append(np.std(peaks))
    noise_sigma.append(options['noise_sigma'])

    stateList = ['I_clamp_GABA_A', 'I_clamp_NMDA']
    maxLag    = 1. / (options['theta_freq'] * 1e-3)
    freq, acval, acVec = extractACStat(stateMonF_e, stateList, maxLag) 

    AC_all.append(acVec)
    freq_mean.append(np.mean(freq))
    freq_std.append(np.std(freq))
    acval_mean.append(np.mean(acval))
    acval_std.append(np.std(acval))

    
###############################################################################
plotSigmas(noise_sigma, peaks_mean, peaks_std, freq_mean, freq_std, acval_mean,
        acval_std)

AC_times = np.arange(len(AC_all[0][0]))
plotACs(AC_all, AC_times, noise_sigma)


###############################################################################
plt.show()
