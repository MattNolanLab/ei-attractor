#!/usr/bin/env python
#
#   plot_EI_velocity.py
#
#   Plot the EI parameter exploration for the bump velocity analysis.
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

from plotting import plot2DTrial
from parameters.param_space import JobTrialSpace2D, DataSpace

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


def aggregate2D(sp, varList):
    varList = ['analysis'] + varList
    return sp.aggregateData(varList, trialNumList='all-at-once', saveData=True)


def computeYX(sp, iterList):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[0][0][0].data['IvelData'], 'net_Ne')
    Ni = DataSpace.getNetParam(sp[0][0][0].data['IvelData'], 'net_Ni')
    return E/Ne, I/Ni


def plotVelSlope(sp, varList, iterList, xlabel="", ylabel="", colorBar=True,
        clBarLabel="", vmin=None, vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    slopes = aggregate2D(sp, varList)
    Y, X = computeYX(sp, iterList)
    return plot2DTrial(X, Y, FR, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)
            
###############################################################################

dirs = \
    ('EI_param_sweep_150pA', (30, 30))

rootDir = "output/velocity/{0}".format(dirs[0])
shape   = dirs[1]

sp = JobTrialSpace2D(shape, rootDir)
iterList  = ['g_AMPA_total', 'g_GABA_total']
    
################################################################################
plt.figure(figsize=(5.1, 2.9))
N = 2
plt.subplot(1, N, 1)

slope = plotVelSlope(sp, ['lineFitSlope'], iterList,
        xlabel="I (nS)",
        ylabel='E (nS)',
        clBarLabel = "Correlation",
        clbarNTicks=3)
###############################################################################
plt.tight_layout()
noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
plt.savefig(sp.rootDir +
    '/../analysis_EI_{0}pA_vel.png'.format(int(noise_sigma)))

