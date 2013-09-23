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
import matplotlib
matplotlib.use('Agg')
from parameters  import JobTrialSpace2D
from figures.EI_plotting import plotACTrial, plotBumpSigmaTrial, \
        plotBumpErrTrial, plotFRTrial
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)



###############################################################################

if (__name__ == '__main__'):
    from matplotlib import rc
    rc('pdf', fonttype=42)
    rc('mathtext', default='regular')
    import matplotlib.pyplot as plt

    # Other
    plt.rcParams['font.size'] = 11
    figSize = (5.1, 3.5)

    dirs = \
        ('150pA',  (31, 31))
        #('300pA',  (31, 31))
        #('0pA',    (31, 31))

    NTrials = 5
    rootDir = "output/no_theta/gamma_bump/{0}".format(dirs[0])
    shape   = dirs[1]

    sp = JobTrialSpace2D(shape, rootDir)
    iterList  = ['g_AMPA_total', 'g_GABA_total']
        
    ################################################################################
    plt.figure(figsize=figSize)
    N = 2
    plt.subplot(1, N, 1)

    acVal = plotACTrial(sp, ['acVal'], iterList,
            trialNumList=range(NTrials),
            xlabel="$g_I$ (nS)",
            ylabel='$g_E$ (nS)',
            clBarLabel = "Correlation",
            clbarNTicks=3,
            vmin=0,
            vmax=0.75,
            r=0, c=1)
    ###############################################################################
    plt.subplot(1, N, 2)
    freq = plotACTrial(sp, ['freq'], iterList,
            trialNumList=range(NTrials),
            xlabel="$g_I$ (nS)",
            clBarLabel = "Frequency (Hz)",
            clbarNTicks=3,
            yticks=False,
            vmin=0,
            vmax=150,
            r=0, c=1)
    ###############################################################################
    plt.tight_layout()
    noise_sigma = sp.getParam(sp[0][1][0].data, 'noise_sigma')
    plt.savefig(sp.rootDir +
            '/../analysis_EI_{0}pA.png'.format(int(noise_sigma)))

    ################################################################################
    # The following is an average over trials
    ################################################################################
    plt.figure(figsize=figSize)
    N = 2
    plt.subplot(1, N, 1)

    bumpSigmaThreshold = 50
    bump_sigma = plotBumpSigmaTrial(sp, ['bump_e', 'sigma'], iterList,
            thr=bumpSigmaThreshold,
            trialNumList=range(NTrials),
            xlabel="$g_I$ (nS)",
            ylabel='$g_E$ (nS)',
            clBarLabel = "Gaussian $\sigma$ (nrns)",
            vmin=0,
            vmax=bumpSigmaThreshold,
            clbarNTicks=4,
            r=0, c=1)
    ###############################################################################
    plt.subplot(1, N, 2)
    bump_err2 = plotBumpErrTrial(sp, ['bump_e', 'err2'], iterList,
            trialNumList=range(NTrials),
            mask=bump_sigma.mask,
            xlabel="$g_I$ (nS)",
            clBarLabel = "Error of fit (Hz)",
            clbarNTicks=3,
            yticks=False,
            vmin=0,
            vmax=200,
            r=0, c=1)
    ###############################################################################
    plt.tight_layout()
    noise_sigma = sp.getParam(sp[0][1][0].data, 'noise_sigma')
    plt.savefig(sp.rootDir +
            '/../analysis_EI_{0}pA_bump.png'.format(int(noise_sigma)))
    ###############################################################################
    ###############################################################################
    plt.figure(figsize=figSize)
    N = 2
    plt.subplot(1, N, 1)

    FR_e = plotFRTrial(sp, ['FR_e', 'avg'], iterList,
            trialNumList=range(NTrials),
            xlabel="$g_I$ (nS)",
            ylabel='$g_E$ (nS)',
            clBarLabel = "E cell firing rate (Hz)",
            clbarNTicks=4,
            vmin=0,
            vmax=50,
            r=0, c=1)
    ###############################################################################
    plt.subplot(1, N, 2)
    FR_threshold=250
    FR_i = plotFRTrial(sp, ['FR_i', 'avg'], iterList,
            trialNumList=range(NTrials),
            thr=FR_threshold,
            xlabel="$g_I$ (nS)",
            clBarLabel = "I cell firing rate (Hz)",
            clbarNTicks=3,
            vmin=0,
            vmax=FR_threshold,
            yticks=False,
            r=0, c=1)
    ###############################################################################
    plt.tight_layout()
    noise_sigma = sp.getParam(sp[0][1][0].data, 'noise_sigma')
    plt.savefig(sp.rootDir +
        '/../analysis_EI_{0}pA_FR.png'.format(int(noise_sigma)))

