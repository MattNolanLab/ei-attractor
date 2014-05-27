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

from EI_plotting          import sweeps, examples
from grid_cell_model.parameters           import JobTrialSpace2D, DataSpace
from grid_cell_model.plotting.global_defs import globalAxesSettings
from EI_plotting.base     import plotOneHist, NoiseDataSpaces
from grid_cell_model.submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--bumpSweep')
parser.add_flag('--slopeSweeps')
parser.add_flag('--slopeErrSweeps')
parser.add_flag('--gammaPowerSweep')
parser.add_flag('--gammaFreqSweep')
parser.add_flag('--gammaExamples')
args = parser.parse_args()


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "panels"
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
gammaRC = [(0, 0), (0, 1), (0, 0)] # (row, col)
gridRC  = [(0, 0), (0, 0), (0, 0)] # (row, col)
bumpDataRoot  = 'output_local/no_theta/gamma_bump'
velDataRoot   = 'output_local/no_theta/velocity'
gridsDataRoot = 'output_local/no_theta/grids'
shape    = (31, 31)



###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

# gamma example rows and columns
exampleRC = ( (5, 5), (15, 15) )

sweepFigSize = (3.5, 2.5)
sweepLeft    = 0.15
sweepBottom  = 0.2
sweepRight   = 0.87
sweepTop     = 0.85
transparent  = True

##############################################################################
gridNTrials = 3
gridTrialNumList = np.arange(gridNTrials)
varList = ['gridnessScore']

grid_cbar_kw= {
    'label'      : 'Gridness score',
    'location'   : 'right',
    'shrink'     : 0.8,
    'pad'        : -0.05,
    'ticks'      : ti.MultipleLocator(0.5),
    'rasterized' : True}

grid_vmin = -0.57
grid_vmax = 1.048

if args.grids or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        sweeps.plotGridTrial(ps.grids[ns_idx], varList, iterList,
                noise_sigma=noise_sigma,
                trialNumList=gridTrialNumList,
                ax=ax,
                r=gridRC[ns_idx][0], c=gridRC[ns_idx][1],
                cbar=cbar, cbar_kw=grid_cbar_kw,
                vmin=grid_vmin, vmax=grid_vmax,
                xlabel='', xticks=False,
                ignoreNaNs=False,
                **kw)
        fname = outputDir + "/suppFigure_no_theta_grids{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()




##############################################################################
sigmaBumpText = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'
sigmaTypes = ['bump_full', 'sigma']
bumpTrialNumList = np.arange(5)
bump_vmin = 0
bump_vmax = 0.516
bump_cbar_kw= dict(
        orientation='vertical',
        shrink=0.8, pad=-0.05,
        ticks=ti.MultipleLocator(0.25),
        label=sigmaBumpText)

if args.bumpSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        sweeps.plotBumpSigmaTrial(ps.bumpGamma[ns_idx], sigmaTypes, iterList,
                noise_sigma=noise_sigma,
                r=gammaRC[ns_idx][0], c=gammaRC[ns_idx][1],
                ax=ax,
                trialNumList=bumpTrialNumList,
                cbar=cbar, cbar_kw=bump_cbar_kw,
                sigmaTitle=False,
                xlabel='', xticks=False,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)

        fname = outputDir + "/suppFigure_no_theta_bump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()


##############################################################################
slopeVarList = ['lineFitSlope']
slope_vmin = 0
slope_vmax = 3.64

slope_cbar_kw= dict(
        location='right',
        shrink = 0.8,
        pad = -0.05,
        labelpad=8,
        label='Slope\n(neurons/s/pA)',
        ticks=ti.MultipleLocator(1))

if args.slopeSweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        _, ax, cax = sweeps.plotVelTrial(ps.v[ns_idx], slopeVarList, iterList,
                noise_sigma, sigmaTitle=False,
                ax=ax,
                cbar=cbar, cbar_kw=slope_cbar_kw,
                vmin=slope_vmin, vmax=slope_vmax,
                xlabel='', xticks=False,
                **kw)
        fname = outputDir + "/suppFigure_no_theta_slope_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
errVarList = ['lineFitErr']
err_vmin = 0
err_vmax = 10

err_cbar_kw = dict(
    label = 'Fit error (neurons/s)',
    orientation = 'vertical',
    shrink = 0.8,
    pad = -0.05,
    ticks = ti.MultipleLocator(4),
    rasterized = True,
    extend='max', extendfrac=0.1)


if args.slopeErrSweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        _, ax, cax = sweeps.plotVelTrial(ps.v[ns_idx], errVarList, iterList,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar=cbar, cbar_kw=err_cbar_kw,
                xlabel='', xticks=False,
                vmin=err_vmin, vmax=err_vmax,
                sigmaTitle=False,
                **kw)
        fname = outputDir + "/suppFigure_no_theta_err_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)



##############################################################################
AC_vmin = -0.21
AC_vmax = 0.96
F_vmin  = 25
F_vmax  = 120

ACVarList = ['acVal']
FVarList  = ['freq']

AC_cbar_kw = dict(
        orientation='vertical',
        ticks=ti.MultipleLocator(0.3),
        shrink=0.8,
        pad=-0.05,
        label='Correlation')
F_cbar_kw = dict(
        orientation='vertical',
        ticks=ti.MultipleLocator(30),
        shrink=0.8,
        pad=-0.05,
        label='Frequency',
        extend='max', extendfrac=0.1)


if args.gammaPowerSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        sweeps.plotACTrial(ps.bumpGamma[ns_idx], ACVarList, iterList,
                noise_sigma=noise_sigma,
                r=gammaRC[ns_idx][0], c=gammaRC[ns_idx][1],
                ax=ax,
                trialNumList=xrange(NTrials),
                xlabel='', xticks=False,
                sigmaTitle=False,
                cbar=cbar, cbar_kw=AC_cbar_kw,
                vmin=AC_vmin, vmax=AC_vmax,
                **kw)
        fname = outputDir + "/suppFigure_no_theta_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()
        

if args.gammaFreqSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        sweeps.plotACTrial(ps.bumpGamma[ns_idx], FVarList, iterList,
                noise_sigma=noise_sigma,
                r=gammaRC[ns_idx][0], c=gammaRC[ns_idx][1],
                ax=ax,
                trialNumList=xrange(NTrials),
                cbar=cbar, cbar_kw=F_cbar_kw,
                sigmaTitle=False,
                vmin=F_vmin, vmax=F_vmax,
                **kw)
        fname = outputDir + "/suppFigure_no_theta_freq_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()
        


##############################################################################
exampleFName = outputDir + "/suppFigure_no_theta_example{0}_{1}.pdf"
exampleTrialNum = 0
exampleFigSize = (2, 0.82)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.82
if args.gammaExamples or args.all:
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
            examples.plotGammaExample(ps.bumpGamma[nsIdx], ax=ax,
                    r=exampleRC[idx][0], c=exampleRC[idx][1],
                    trialNum=exampleTrialNum,
                    tStart = 2e3, tEnd=2.25e3,
                    noise_sigma=nsAnn)
            plt.savefig(fname, dpi=300, transparent=True)
            plt.close()



