#!/usr/bin/env python
#
#   suppFigure_grids_vs_bumps.py
#
#   Supplementary figure: scatter plots of gridness score vs. bump width.
#
import math
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.transforms import Bbox
from copy import deepcopy

import default_settings as ds
from EI_plotting          import sweeps, scatter
from EI_plotting          import aggregate as aggr
from EI_plotting.base     import NoiseDataSpaces
from grid_cell_model.parameters           import JobTrialSpace2D
from grid_cell_model.plotting.global_defs import globalAxesSettings
from submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--scatterPlot')
parser.add_argument('--expScatter', action='store_true')
args = parser.parse_args()

outputDir = "output_figures"
ps = ds.getDefaultParamSpaces()


##############################################################################
# Scatter plot of gridness score vs. bump sigma^-1
scatterFigSize = (8.27, 11.69)
scatterLeft   = 0.12
scatterBottom = 0.17
scatterRight  = 0.98
scatterTop    = 0.92
scatterTransparent = True

scatterColorFigSize = (1.5, 1.5)

ignoreNaNs = True


class CustomScale(mscale.ScaleBase):
    name = 'custom'

    def __init__(self, axis, **kwargs):
        mscale.ScaleBase.__init__(self)
        self.thresh = None #thresh

    def get_transform(self):
        return self.CustomTransform(self.thresh)

    def set_default_locators_and_formatters(self, axis):
        pass

    class CustomTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            return 10**(a*2)

        def inverted(self):
            return CustomScale.InvertedCustomTransform(self.thresh)

    class InvertedCustomTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            return log10(a)

        def inverted(self):
            return CustomScale.CustomTransform(self.thresh)


mscale.register_scale(CustomScale)

##############################################################################
xlabel = 'P(bumps)'
ylabel = 'Gridness score'

if args.scatterPlot or args.all:
    fig = plt.figure(figsize=scatterFigSize)
    ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
        scatterTop))

    isBumpData = []
    gridData = []
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        isBumpData.append(aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList,
            ignoreNaNs=True))
        gridData.append(aggr.GridnessScore(ps.grids[ns_idx], ds.iterList,
            ignoreNaNs=True))

    fig = plt.figure(figsize=scatterFigSize)
    scatterPlot = scatter.FullScatterPlot(
            isBumpData, gridData, None, None, ds.iterList, None, None,
            s=25,
            linewidth=0.3,
            color2D=True,
            xlabel=xlabel,
            ylabel=ylabel,
            sigmaTitle=True,
            noise_sigmas=ps.noise_sigmas,
            ignoreNaNs=ignoreNaNs,
            captionLetters=('A', 'B', 'C'),
            fig=fig)
    scatterPlot.plot(captionLeft=-0.1, plotcolorbar=False)
    l = 0.14
    w = 0.165
    scatterPlot.plotColorbar(left=l, bottom=.85, right=l+w, top=.95)
    scatterPlot.set_titleSizes(16)
    if args.expScatter:
        for ns_idx, _ in enumerate(ps.noise_sigmas):
            ax = scatterPlot.axes[ns_idx]
            ax.set_xscale('custom')
            ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
            ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
            ax.set_xlim([-0.3, 1.002])
        fname = outputDir + "/suppFigure_grids_vs_bumps_exp.pdf"
    else:
        fname = outputDir + "/suppFigure_grids_vs_bumps.pdf"

    fig.savefig(fname, dpi=300)

