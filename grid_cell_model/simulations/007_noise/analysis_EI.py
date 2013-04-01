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
import matplotlib.pyplot as plt

from analysis.visitors    import AutoCorrelationVisitor
from parameters           import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings, createColorbar
import logging as lg
lg.basicConfig(level=lg.DEBUG)


# Other
plt.rcParams['font.size'] = 12

###############################################################################

def plot2DTrial(sp, xlabel="", ylabel="", colorBar=True, title="",
        clbarNTicks=2):
    shape = sp.getShape()
    rows = shape[0]
    cols = shape[1]
    trialNum = 0
    acVal = np.ndarray(shape)
    for r in xrange(rows):
        for c in xrange(cols):
            if (len(sp[r][c]) == 0):
                acVal[r][c] = 0
            else:
                data = sp[r][c][trialNum].data['analysis']
                acVal[r][c] = np.mean(data['acVal'])

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.pcolormesh(acVal)
    createColorbar(ax, acVal, 'Correlation', nticks=clbarNTicks)
    if (xlabel != ""):
        plt.xlabel(xlabel)
    if (ylabel != ""):
        plt.ylabel(ylabel)

    return acVal


###############################################################################

#rootDir = 'output_local/2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample'
#shape = (2, 2)

rootDir = 'output_local/2013-03-30T19-34-44_EI_param_sweep_0pA_full'
shape = (20, 20)
sp = JobTrialSpace2D(shape, rootDir)


monName   = 'stateMonF_e'
stateList = ['I_clamp_GABA_A']
forceUpdate = False
visitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=forceUpdate)
sp.visit(visitor)

###############################################################################

acVal = plot2DTrial(sp,
        xlabel="E coupling strength (nS) TODO",
        ylabel='I coupling strength (nS) TODO',
        clbarNTicks=5)

plt.show()
