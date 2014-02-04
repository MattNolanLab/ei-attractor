#!/usr/bin/env python
#
#   figure1.py
#
#   Noise publication Figure 1.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from EI_plotting      import sweeps, examples, details, segmentation
from EI_plotting      import aggregate as aggr
from EI_plotting.base import NoiseDataSpaces, getOption, plotStateSignal
from parameters       import JobTrialSpace2D
from data_storage     import DataStorage
from data_storage.sim_models.ei import extractSummedSignals
import plotting.low_level
from plotting.global_defs import prepareLims
from analysis         import clustering
import flagparse

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=3
gridTrialNumList = np.arange(NTrials)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas   = [0, 150, 300]
exampleIdx     = [(1, 22), (1, 22), (1, 22)] # (row, col)
bumpDataRoot   = None
velDataRoot    = None
gridsDataRoot  = 'output_local/even_spacing/grids'
singleDataRoot = 'output_local/single_neuron'
shape = (31, 31)

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--examplesFlag')
parser.add_flag('--detailed_noise')
parser.add_flag('--Vm_examples')
args = parser.parse_args()

##############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True

cbar_kw= {
    'label'      : 'Gridness score',
    'location'   : 'right',
    'shrink'     : 0.8,
    'pad'        : -0.05,
    'ticks'      : ti.MultipleLocator(0.5),
    'rasterized' : True}

vmin = -1.1
vmax = 1.1

##############################################################################
exampleRC = ( (5, 15), (15, 5) )
sliceAnn = None
grids_cmap = 'RdBu_r'

ann0 = dict(
        txt='b',
        rc=exampleRC[0],
        xytext_offset=(1.5, 1),
        color='black')
ann1 = dict(
        txt='a',
        rc=exampleRC[1],
        xytext_offset=(0.5, 1.5),
        color='black')
ann = [ann0, ann1]


varList = ['gridnessScore']

if args.grids or args.all:
    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotGridTrial(ps.grids[0], varList, iterList, ps.noise_sigmas[0],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[0][0], c=exampleIdx[0][1],
            cbar=False, cbar_kw=cbar_kw,
            cmap=grids_cmap,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            annotations=ann,
            sliceAnn=sliceAnn)
    fname = outputDir + "/grids_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


    # noise_sigma = 150 pA
    #for a in sliceAnn:
    #    a['letterColor'] = 'black'
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotGridTrial(ps.grids[1], varList, iterList, ps.noise_sigmas[1],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[1][0], c=exampleIdx[1][1],
            cbar=False, cbar_kw=cbar_kw,
            cmap=grids_cmap,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            ylabel='', yticks=False,
            annotations=ann,
            sliceAnn=sliceAnn)
    fname = outputDir + "/grids_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    _, _, cax = sweeps.plotGridTrial(ps.grids[2], varList, iterList, ps.noise_sigmas[2],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[2][0], c=exampleIdx[2][1],
            cbar_kw=cbar_kw,
            cmap=grids_cmap,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            ylabel='', yticks=False,
            annotations=ann,
            sliceAnn=sliceAnn)
    #for label in cax.yaxis.get_ticklabels():
    #    label.set_ha('right')
    #cax.tick_params(pad=30)
    fname = outputDir + "/grids_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


##############################################################################
# Grid field examples
exampleGridFName = outputDir + "/grids_examples_{0}pA_{1}.pdf"
exampleACFName = outputDir + "/grids_examples_acorr_{0}pA_{1}.pdf"
exTransparent = True
exampleFigSize = (1, 1.2)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.85

if args.examplesFlag or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            # Grid field
            fname = exampleGridFName.format(noise_sigma, idx)
            fig = plt.figure(figsize=exampleFigSize)
            gs = examples.plotOneGridExample(ps.grids[ns_idx], rc, iterList,
                    exIdx=exampleIdx[idx],
                    xlabel=False, ylabel=False,
                    xlabel2=False, ylabel2=False, 
                    maxRate=True, plotGScore=False,
                    fig=fig)
            gs.update(left=exampleLeft, bottom=exampleBottom, right=exampleRight,
                    top=exampleTop)
            plt.savefig(fname, dpi=300, transparent=exTransparent)
            plt.close()

            # Autocorrelation
            fname = exampleACFName.format(noise_sigma, idx)
            fig= plt.figure(figsize=exampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                exampleRight, exampleTop))
            gs = examples.plotOneGridACorrExample(ps.grids[ns_idx], rc, ax=ax)
            plt.savefig(fname, dpi=300, transparent=exTransparent)
            plt.close()
    

##############################################################################
# Membrane potential examples

def openJob(rootDir, noise_sigma):
    fileTemplate = "noise_sigma{0}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(int(noise_sigma))
    return DataStorage.open(fileName, 'r')

def drawVm(data, noise_sigma, xScaleBar=None, yScaleBar=None,
        ax=plt.gca(), sigmaTitle=True):
    yScaleX = 0.5
    yScaleY = 1.1
    yScaleXOffset = 0.06
    scaleTextSize = 'x-small'

    plotTStart = 5e3
    plotTEnd   = 5.25e3

    stateYlim = [-80, -40]

    theta_start_t = getOption(data, 'theta_start_t')
    #theta_start_t = 1e3
    simTime = getOption(data, 'time')

    mon_e = data['stateMon_e']

    # E cell Vm
    t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax, t, VmMiddle, labely='', color='red',
            scaleBar=xScaleBar, scaleText="ms", scaleX=0.5, scaleY=1.1,
            scaleTextSize=scaleTextSize)
    if (yScaleBar is not None):
        plotting.low_level.yScaleBar(yScaleBar, yScaleX, yScaleY,
                ax=ax,
                unitsText='mV', textXOffset=yScaleXOffset,
                size='x-small')
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.ylim(stateYlim)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), loc='left',
                size=scaleTextSize, y=0.9)


VmExampleFigSize = (2.5, 1.25)
VmExampleLeft   = 0.01
VmExampleBottom = 0.01
VmExampleRight  = 0.999
VmExampleTop    = 0.6
VmExampleXScalebar = 50 # ms
VmExampleYScalebar = 10 # mV

if args.Vm_examples or args.all:
    for ns_idx, noise_sigma in enumerate(noise_sigmas):
        fig = plt.figure(figsize=VmExampleFigSize)
        ax = fig.add_axes(Bbox.from_extents(VmExampleLeft, VmExampleBottom,
            VmExampleRight, VmExampleTop))
        ds = openJob(singleDataRoot, noise_sigma)
        kw = {}
        if (ns_idx == 2):
            kw['xScaleBar'] = VmExampleXScalebar
            kw['yScaleBar'] = VmExampleYScalebar
        drawVm(ds, noise_sigma=noise_sigma, ax=ax, **kw)

        fname = outputDir + "/grids_Vm_example_{0}.pdf"
        plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Detailed noise plots
EI13Root  = 'output_local/detailed_noise/grids/EI-1_3'
EI31Root  = 'output_local/detailed_noise/grids/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 1


detailFigSize = (4.5, 2.8)
detailLeft   = 0.15
detailBottom = 0.2
detailRight  = 0.95
detailTop    = 0.8
if args.detailed_noise or args.all:
    ylabelPos = -0.15

    types = ('grids', 'gridnessScore')
    fig = plt.figure(figsize=detailFigSize)
    ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
        detailTop))
    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            color='red', markerfacecolor='red', zorder=10)
    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            ylabel='Gridness score', ylabelPos=ylabelPos,
            color='#505050')
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))
    ax.set_ylim([-0.6, 1.2])
    leg = ['a',  'b']
    l = ax.legend([p31, p13], leg, loc=(0.8, 1), fontsize='small', frameon=False,
            numpoints=1, handletextpad=0.05)

    fname = outputDir + "/grids_detailed_noise_gscore.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



