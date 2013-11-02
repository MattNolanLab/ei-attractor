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
import os
import numpy as np
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, xlabel, ylabel, hold, hist, \
        subplot
from matplotlib.ticker  import MaxNLocator, AutoMinorLocator
import numpy.ma          as ma

from plotting.global_defs import globalAxesSettings
from figures.EI_plotting import plot2DTrial
from parameters.param_space import JobTrialSpace2D, DataSpace
from figures.figures_shared import plotOneHist
from figures.EI_plotting    import aggregate2D, plotVelTrial

import logging as lg
lg.basicConfig(level=lg.WARN)
#lg.basicConfig(level=lg.INFO)

dir = ('{0}pA', (31, 31))


def getAll(sp, varList, iterList):
    vlSlope = varList + ['lineFitSlope']
    vlErrs = varList + ['lineFitErr']
    slopes = np.abs(aggregate2D(sp, vlSlope, funReduce=np.sum))
    errs   = np.abs(aggregate2D(sp, vlErrs, funReduce=np.sum))
    Y, X = computeYX(sp, iterList)
    return ((X, Y), slopes.flatten(), errs.flatten())

def plotSlopeVsError(sp, varList, iterList):
    ((X, Y), slopes, errs) = getAll(sp, varList, iterList)

    ax = plt.gca()
    globalAxesSettings(ax)
    plot(errs, slopes, 'o', markersize=3)
    xlabel("Error of fit (nrns/s)")
    ylabel("Slope (nrnrs/s/pA)")
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))


def plotHistograms(sp, varList, iterList):
    ((X, Y), slopes, errs) = getAll(sp, varList, iterList)
    slopes[np.isnan(slopes)] = []
    errs[np.isnan(errs)] = []
    fig = plt.gcf()
    ax = subplot(1, 2, 1)
    plotOneHist(slopes)
    xlabel("Slope (nrnrs/s/pA)")
    ylabel("Count")

    ax = subplot(1, 2, 2)
    plotOneHist(errs)
    xlabel("Error of fit (nrns/s)")
    




            
###############################################################################
noise_sigmas = [0, 150, 300]


for noise_sigma in noise_sigmas:
    file = dir[0].format(int(noise_sigma))
    rootDir = "output/even_spacing/velocity/{0}".format(file)
    shape   = dir[1]

    sp = JobTrialSpace2D(shape, rootDir)
    iterList  = ['g_AMPA_total', 'g_GABA_total']

        
    ################################################################################
    plt.figure(figsize=(5.1, 2.9))
    N = 2
    plt.subplot(1, N, 1)
    slope = plotVelTrial(sp, ['lineFitSlope'], iterList,
            xlabel="I (nS)",
            ylabel='E (nS)',
            clBarLabel = "Slope",
            clbarNTicks=3)

    plt.subplot(1, N, 2)
    lineFitErr = plotVelTrial(sp, ['lineFitErr'], iterList,
            xlabel="I (nS)",
            ylabel='E (nS)',
            clBarLabel = "Line fit error",
            clbarNTicks=4,
            vmin=0,
            vmax=18)

    plt.tight_layout()
    noise_sigma = sp.getParam(sp[0][0][0].data['IvelData'][0], 'noise_sigma')
    plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA_vel.png'.format(int(noise_sigma)))
    ###############################################################################
    plt.figure(figsize=(3, 3))
    plotSlopeVsError(sp, [], iterList)
    plt.tight_layout()
    plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA_slope_vs_err.png'.format(int(noise_sigma)))
    ###############################################################################
    plt.figure(figsize=(5.1, 2.9))
    plotHistograms(sp, [], iterList)
    plt.tight_layout()
    plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA_vel_hist.png'.format(int(noise_sigma)))
