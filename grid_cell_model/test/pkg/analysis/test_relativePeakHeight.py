#! /usr/bin/env python
#
#   test_relativePeakHeight.py
#
#   Test a function.
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
from analysis.signal import localExtrema, relativePeakHeight, autoCorrelation,\
        butterBandPass

from matplotlib.pyplot import *


def readTxtFile(fileName):
    res = []
    f = open(dFile, 'r')
    for line in f:
        res.append(float(line))
    return np.array(res)


def processPeaks(sig, plotPeaks = False, plotHeights=False):
    ext_idx, ext_t = localExtrema(sig)
    ext = sig[ext_idx]
    times = np.arange(sig.size)
    heights = relativePeakHeight(ext, np.min)

    if (plotPeaks):
        figure()
        plot(times, sig)
        hold('on')
        plot(times[ext_idx], sig[ext_idx], '.')
        ylim([-1, 1])

    if (plotHeights):
        # plot the heights as lines
        for h_idx in range(len(heights)):
            h = heights[h_idx]
            d = sig[ext_idx[h_idx]]
            plot([times[ext_idx[h_idx]]] * 2, [d, d - h*ext_t[h_idx]])

    return (ext_idx, ext_t, heights)


dFile = 'auto_correlation.txt'
data = readTxtFile(dFile)


ext_idx, ext_t, heights = processPeaks(data, plotPeaks=True, plotHeights=True)
ylabel('Normalized autocorrelation')
xlabel('Time (ms)')
title('Peak finding and height detection')

show()
