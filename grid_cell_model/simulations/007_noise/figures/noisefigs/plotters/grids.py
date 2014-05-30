from __future__ import absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from grid_cell_model.parameters       import JobTrialSpace2D
from grid_cell_model.data_storage     import DataStorage
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals
import grid_cell_model.plotting.low_level as low_level

from ..EI_plotting import sweeps, examples, details, scatter
from ..EI_plotting import aggregate as aggr
from ..EI_plotting.base import getOption, plotStateSignal
from ..EI_plotting import scaling
from .base import FigurePlotter, SweepPlotter

__all__ = [
    'GridSweepsPlotter',
    'GridExamplesPlotter',
    'VMExamplesPlotter',
    'GridDetailedNoisePlotter',
    'GridsDiffSweep',
    'GridBumpScatterPlotter',
]
         

##############################################################################

class GridSweepsPlotter(SweepPlotter):
    cmap = 'jet'
    varList = ['gridnessScore']

    def __init__(self, *args, **kwargs):
        super(GridSweepsPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        example_idx = self.config['grids']['example_idx']
        trial_num_list = np.arange(self.config['grids']['ntrials'])

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/grids_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx == 2:
                    kw['cbar'] = True
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                sweeps.plotGridTrial(
                        ps.grids[ns_idx],
                        self.varList,
                        self.config['iter_list'],
                        ps.noise_sigmas[ns_idx],
                        trialNumList=trial_num_list,
                        r=example_idx[ns_idx][0], c=example_idx[ns_idx][1],
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        cmap=self.cmap,
                        vmin=myc['vmin'], vmax=myc['vmax'],
                        ignoreNaNs=True,
                        annotations=myc['ann'],
                        sliceAnn=None,
                        **kw)


###############################################################################
## Grid field examples
class GridExamplesPlotter(FigurePlotter):
    exampleGridFName = "/grids_examples_{0}pA_{1}.pdf"
    exampleACFName = "/grids_examples_acorr_{0}pA_{1}.pdf"

    def __init__(self, *args, **kwargs):
        super(GridExamplesPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        ps = self.env.ps
        exampleLeft, exampleBottom, exampleRight, exampleTop = myc['ax_box']
        example_idx = self.config['grids']['example_idx']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            for idx, rc in enumerate(self.config['grids']['example_rc']):
                # Grid field
                fname = (self.config['output_dir'] + 
                            self.exampleGridFName.format(noise_sigma, idx))
                fig = self._get_final_fig(myc['fig_size'])
                gs = examples.plotOneGridExample(
                        ps.grids[ns_idx],
                        rc,
                        self.config['iter_list'],
                        exIdx=example_idx[idx],
                        xlabel=False, ylabel=False,
                        xlabel2=False, ylabel2=False, 
                        maxRate=True, plotGScore=False,
                        fig=fig)
                gs.update(left=exampleLeft, bottom=exampleBottom, right=exampleRight,
                        top=exampleTop)
                plt.savefig(fname, dpi=300, transparent=myc['transparent'])
                plt.close()
    
                # Autocorrelation
                fname = (self.config['output_dir'] +
                            self.exampleACFName.format(noise_sigma, idx))
                fig= plt.figure(figsize=myc['fig_size'])
                ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                    exampleRight, exampleTop))
                gs = examples.plotOneGridACorrExample(ps.grids[ns_idx], rc, ax=ax)
                plt.savefig(fname, dpi=300, transparent=myc['transparent'])
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
        low_level.yScaleBar(yScaleBar, yScaleX, yScaleY,
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

class VMExamplesPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VMExamplesPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(VmExampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(VmExampleLeft, VmExampleBottom,
                VmExampleRight, VmExampleTop))
            ds = openJob(self.config['singleDataRoot'], noise_sigma)
            kw = {}
            if (ns_idx == 2):
                kw['xScaleBar'] = VmExampleXScalebar
                kw['yScaleBar'] = VmExampleYScalebar
            drawVm(ds, noise_sigma=noise_sigma, ax=ax, **kw)

            fname = self.config['output_dir'] + "/grids_Vm_example_{0}.pdf"
            plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
            plt.close()


##############################################################################
# Detailed noise plots
EI13Root  = 'simulation_data/submission/detailed_noise/grids/EI-1_3'
EI31Root  = 'simulation_data/submission/detailed_noise/grids/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 1


detailFigSize = (4.5, 2.8)
detailLeft   = 0.15
detailBottom = 0.2
detailRight  = 0.95
detailTop    = 0.8

class GridDetailedNoisePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GridDetailedNoisePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ylabelPos = -0.15

        types = ('grids', 'gridnessScore')
        fig = self._get_final_fig(detailFigSize)
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

        fname = self.config['output_dir'] + "/grids_detailed_noise_gscore.pdf"
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()


##############################################################################
# Parameter sweep of the difference between noise_150 and noise_0
gridTypes = ['grids', 'gridnessScore']

class GridsDiffSweep(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(GridsDiffSweep, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        example_idx = self.config['grids']['example_idx']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
            fname = (self.config['output_dir'] +
                     "/grids_diff_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                which = ns_idx
                sweeps.plotDiffTrial(
                    ps.grids, self.config['iter_list'], which,
                    self.config['grids']['ntrials'], gridTypes,
                    ax=ax,
                    ignoreNaNs=True,
                    r=example_idx[0][0], c=example_idx[0][1],
                    cbar=True, cbar_kw=myc['cbar_kw'],
                    symmetricLimits=True,
                    cmap='RdBu_r'
                )


##############################################################################
# Scatter plot of gridness score vs. P(bump)
##############################################################################
class GridBumpScatterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GridBumpScatterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        output_dir = self.config['output_dir']

        scatterFigSize = (8.27, 11.69)
        scatterLeft   = 0.12
        scatterBottom = 0.17
        scatterRight  = 0.98
        scatterTop    = 0.92
        scatterTransparent = True
        
        scatterColorFigSize = (1.5, 1.5)
        
        ignoreNaNs = True

        xlabel = 'P(bumps)'
        ylabel = 'Gridness score'

        fig = self._get_final_fig(scatterFigSize)
        ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
            scatterTop))

        isBumpData = []
        gridData = []
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            isBumpData.append(aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                ignoreNaNs=True, normalizeTicks=False))
            gridData.append(aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                ignoreNaNs=True, normalizeTicks=False))

        fig = self._get_final_fig(scatterFigSize)
        scatterPlot = scatter.FullScatterPlot(
                isBumpData, gridData, None, None, iter_list, None, None,
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

        # Normal scale
        fname = output_dir + "/suppFigure_grids_vs_bumps.pdf"
        fig.savefig(fname, dpi=300)
        
        # Exponential scale
        for ns_idx, _ in enumerate(ps.noise_sigmas):
            ax = scatterPlot.axes[ns_idx]
            ax.set_xscale('exponential')
            ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
            ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
            ax.set_xlim([-0.3, 1.002])
        fname = output_dir + "/suppFigure_grids_vs_bumps_exp.pdf"
        fig.savefig(fname, dpi=300)

