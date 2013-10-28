#!/usr/bin/env python
#
#   suppFigure_no_theta.py
#
#   Supplementary figure showing the results without theta stimulation.
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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy


import EI_plotting as EI
from parameters  import JobTrialSpace2D, DataSpace
from plotting.global_defs import globalAxesSettings
from figures_shared import plotOneHist, NoiseDataSpaces
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "."
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
gammaRC    = [(0, 0), (0, 1), (0, 0)] # (row, col)
bumpDataRoot  = 'output_local/no_theta/gamma_bump'
velDataRoot   = None
gridsDataRoot = None
shape    = (31, 31)

gammaSweep = 1
examples   = 1


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

# gamma example rows and columns
exampleRC = ( (5, 5), (15, 15) )


sweepFigSize = (2, 2.8)
sweepLeft    = 0.17
sweepBottom  = 0.1
sweepRight   = 0.95
sweepTop     = 0.9
transparent  = True

AC_vmin = -0.21
AC_vmax = 0.96
F_vmin  = 25
F_vmax  = 120

ACVarList = ['acVal']
FVarList  = ['freq']

AC_cbar_kw = dict(
        orientation='horizontal',
        ticks=ti.MultipleLocator(0.3),
        shrink=0.8,
        pad=0.2,
        label='Correlation')
F_cbar_kw = dict(
        orientation='horizontal',
        ticks=ti.MultipleLocator(30),
        shrink=0.8,
        pad=0.2,
        label='Frequency',
        extend='max', extendfrac=0.1)


ann_color = 'black'
ann0 = dict(
        txt='B',
        rc=exampleRC[0],
        xytext_offset=(1.5, 0),
        color=ann_color)
ann1 = dict(
        txt='C',
        rc=exampleRC[1],
        xytext_offset=(1.5, 1),
        color=ann_color)
ann = [ann0, ann1]
annF = [deepcopy(ann0), deepcopy(ann1)]
annF[0]['txt'] = 'B'


if (gammaSweep):
    # noise_sigma = 0 pA
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[0], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            r=gammaRC[0][0], c=gammaRC[0][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            xlabel='', xticks=False,
            cbar=False, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/suppFigure_no_theta_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[0], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            r=gammaRC[0][0], c=gammaRC[0][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            xlabel='', xticks=False,
            ylabel='', yticks=False,
            cbar=False, cbar_kw=F_cbar_kw,
            sigmaTitle=True,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/suppFigure_no_theta_freq_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        

    # noise_sigma = 150 pA
    for a in ann:
        a['color'] = 'black'
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[1], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[1],
            r=gammaRC[1][0], c=gammaRC[1][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            xlabel='', xticks=False,
            cbar=False, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/suppFigure_no_theta_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[1], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[1],
            r=gammaRC[1][0], c=gammaRC[1][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            xlabel='', xticks=False,
            ylabel='', yticks=False,
            cbar=False, cbar_kw=F_cbar_kw,
            sigmaTitle=True,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/suppFigure_no_theta_freq_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        

    # noise_sigma = 300 pA
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[2], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[2],
            r=gammaRC[2][0], c=gammaRC[2][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            cbar=True, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/suppFigure_no_theta_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotACTrial(ps.bumpGamma[2], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[2],
            r=gammaRC[2][0], c=gammaRC[2][1],
            ax=ax,
            trialNumList=xrange(NTrials),
            ylabel='', yticks=False,
            cbar=True, cbar_kw=F_cbar_kw,
            sigmaTitle=True,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/suppFigure_no_theta_freq_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        


##############################################################################
exampleFName = outputDir + "/suppFigure_no_theta_example{0}_{1}.pdf"
exampleTrialNum = 0
exampleFigSize = (2, 0.82)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.82
if (examples):
    for nsIdx, ns in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            fname = exampleFName.format(ns, idx)
            fig = plt.figure(figsize=exampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                exampleRight, exampleTop))
            if (idx == 0):
                nsAnn = ns
            else:
                nsAnn = None
            EI.plotGammaExample(ps.bumpGamma[nsIdx], ax=ax,
                    r=exampleRC[idx][0], c=exampleRC[idx][1],
                    trialNum=exampleTrialNum,
                    tStart = 2e3, tEnd=2.25e3,
                    noise_sigma=nsAnn)
            plt.savefig(fname, dpi=300, transparent=True)
            plt.close()



