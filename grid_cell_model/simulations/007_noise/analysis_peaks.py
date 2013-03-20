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

from data_storage           import DataStorage
from analysis.signal        import splitSigToThetaCycles, globalExtrema
from plotting.global_defs   import globalAxesSettings


# Other
figSize = (12,8)
plt.rcParams['font.size'] = 16


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
        peaks.append(globalExtrema(sig_split, extFunc))

    return np.array(peaks)



###############################################################################


output_dir = 'output_local'
fileNamePrefix = ''

peaks_mean = []
peaks_std  = []
noise_sigma = []

for jobNum in jobNums:
    in_fname = "{0}/{1}job{2:05}_output".format(output_dir, fileNamePrefix, jobNum)
    output_fname = in_fname
    d = DataStorage.open(in_fname + ".h5", 'r')
    options = d['options']
    ei_net  = d['ei_net']
    
    
    peaks = peakAmp(
            stateMon = d['stateMonF_e'],
            startT   = options['theta_start_t'],
            thetaFreq = options['theta_freq'],
            fieldStr = 'I_clamp_GABA_A',
            extFunc = np.max).flat
    peaks_mean.append(np.mean(peaks))
    peaks_std.append(np.std(peaks))
    noise_sigma.append(options['noise_sigma'])
    
plt.figure()
plt.errorbar(noise_sigma, peaks_mean, peaks_std, fmt='o-')
globalAxesSettings(plt.gca())
plt.margins(0.05)
plt.gca().xaxis.set_ticks(noise_sigma)
plt.ion()
plt.show()

