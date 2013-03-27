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

from data_storage           import DataStorage
from analysis.signal        import splitSigToThetaCycles, globalExtremum, corr,\
        localExtrema, butterBandPass, autoCorrelation, downSample
from plotting.global_defs   import globalAxesSettings

import matplotlib
print matplotlib.__version__, matplotlib.__path__


# Other
figSize = (7,5)
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

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


###############################################################################


output_dir = 'output_local'
fileNamePrefix = ''

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
    d = DataStorage.open(in_fname + ".h5", 'r')
    options     = d['options']
    ei_net      = d['ei_net']
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

    freq   = [] # Frequency of input signal
    acval  = [] # Auto-correlation at the corresponding frequency
    for n_id in range(len(stateMonF_e)):
    #for n_id in range(5):
        print "n_id: ", n_id
        sig = stateMonF_e[n_id]['events']['I_clamp_GABA_A']
        times = stateMonF_e[n_id]['events']['times']
        sig, times = downSample(sig, 10, times)
        dt = times[1] - times[0]
        sig = butterBandPass(sig, dt*1e-3, 20, 200)
        slice = 2. / (options['theta_freq'] * 1e-3) / dt
        ac = autoCorrelation(sig - np.mean(sig), max_lag=slice-1, norm=True)
        ext_idx, ext_t = localExtrema(ac)

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

    freq_mean.append(np.mean(freq))
    freq_std.append(np.std(freq))
    acval_mean.append(np.mean(acval))
    acval_std.append(np.std(acval))

    
plt.figure(figsize=figSize)

plt.subplot2grid((2, 2), (0, 0))
plotNoiseSigma(noise_sigma, peaks_mean, peaks_std,
        xlabel = "",
        ylabel="E cell max. $I_{syn}$ / $\\theta$ \n (pA)",
        title="A",
        newFigure=False)
plt.gca().set_xticklabels([])

plt.subplot2grid((2, 2), (1, 0))
plotNoiseSigma(noise_sigma, freq_mean, freq_std,
        ylabel="Frequency (Hz)",
        title="B",
        newFigure=False)

plt.subplot2grid((2, 2), (1, 1))
plotNoiseSigma(noise_sigma, acval_mean, acval_std, 
        ylabel="Mean\nauto-correlation",
        title="C",
        newFigure=False)
plt.tight_layout(w_pad=0.0, h_pad=2.5)

plt.savefig(output_dir + "/peak_analysis.pdf")

plt.show()
