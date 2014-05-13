#!/usr/bin/env python
'''
Supplementary figure with parameter sweeps of grid fields when velocity input
is switched OFF, but place cell input is left ON.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from EI_plotting      import sweeps, examples, details, segmentation
from EI_plotting      import aggregate as aggr
from EI_plotting.base import NoiseDataSpaces, getOption, plotStateSignal
from grid_cell_model.parameters       import JobTrialSpace2D
from grid_cell_model.data_storage     import DataStorage
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals
import plotting.low_level
from grid_cell_model.plotting.global_defs import prepareLims
from grid_cell_model.analysis import clustering
from submitting import flagparse

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
gridsDataRoot  = 'output_local/even_spacing/grids_no_velocity'
shape = (31, 31)

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--examplesFlag')
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
    'ticks'      : ti.MultipleLocator(0.3),
    'rasterized' : True}

vmax = 0.622
vmin = -0.622

##############################################################################
exampleRC = ( (5, 15), (15, 5) )
sliceAnn = None
grids_cmap = 'jet'

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
    ns_idx = 1
    noise_sigma = 150
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotGridTrial(ps.grids[ns_idx], varList, iterList, noise_sigma,
            trialNumList=gridTrialNumList,
            ax=ax,
            cbar=True, cbar_kw=cbar_kw,
            cmap=grids_cmap,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            annotations=ann,
            sliceAnn=sliceAnn)
    fname = outputDir + "/grids_no_velocity_sweeps{0}.pdf"
    fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
    plt.close()



###############################################################################
## Grid field examples
#exampleGridFName = outputDir + "/grids_examples_{0}pA_{1}.pdf"
#exampleACFName = outputDir + "/grids_examples_acorr_{0}pA_{1}.pdf"
#exTransparent = True
#exampleFigSize = (1, 1.2)
#exampleLeft   = 0.01
#exampleBottom = 0.01
#exampleRight  = 0.99
#exampleTop    = 0.85
#
#if args.examplesFlag or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        for idx, rc in enumerate(exampleRC):
#            # Grid field
#            fname = exampleGridFName.format(noise_sigma, idx)
#            fig = plt.figure(figsize=exampleFigSize)
#            gs = examples.plotOneGridExample(ps.grids[ns_idx], rc, iterList,
#                    exIdx=exampleIdx[idx],
#                    xlabel=False, ylabel=False,
#                    xlabel2=False, ylabel2=False, 
#                    maxRate=True, plotGScore=False,
#                    fig=fig)
#            gs.update(left=exampleLeft, bottom=exampleBottom, right=exampleRight,
#                    top=exampleTop)
#            plt.savefig(fname, dpi=300, transparent=exTransparent)
#            plt.close()
#
#            # Autocorrelation
#            fname = exampleACFName.format(noise_sigma, idx)
#            fig= plt.figure(figsize=exampleFigSize)
#            ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
#                exampleRight, exampleTop))
#            gs = examples.plotOneGridACorrExample(ps.grids[ns_idx], rc, ax=ax)
#            plt.savefig(fname, dpi=300, transparent=exTransparent)
#            plt.close()
#    
#
