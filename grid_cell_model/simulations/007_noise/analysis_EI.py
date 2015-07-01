#!/usr/bin/env python
#
'''
Perform analysis on whole 2D data sets.
'''
import time
import matplotlib
matplotlib.use('agg')

import grid_cell_model.visitors as vis
import grid_cell_model.visitors.spikes
import grid_cell_model.visitors.bumps
import grid_cell_model.visitors.signals
import grid_cell_model.visitors.plotting
import grid_cell_model.visitors.plotting.spikes
import grid_cell_model.visitors.plotting.grids
import grid_cell_model.visitors.plotting.grids_ipc
from grid_cell_model.parameters import JobTrialSpace2D
from grid_cell_model.submitting import flagparse

import common.analysis as common

###############################################################################
parser = flagparse.FlagParser()
parser.add_argument('--row',          type=int, required=True)
parser.add_argument('--col',          type=int, required=True)
parser.add_argument('--shapeRows',    type=int, required=True)
parser.add_argument('--shapeCols',    type=int, required=True)
parser.add_argument('--forceUpdate',  type=int, required=True)
parser.add_argument("--output_dir",   type=str, required=True)
parser.add_argument("--job_num",      type=int) # unused
parser.add_argument("--type",         type=str, choices=common.allowedTypes, required=True, nargs="+")
parser.add_argument("--bumpSpeedMax", type=float)

o = parser.parse_args()

###############################################################################
startT = time.time()

shape = (o.shapeRows, o.shapeCols)
dataPoints = [(o.row, o.col)]
trialNums = None

sp = JobTrialSpace2D(shape, o.output_dir, dataPoints=dataPoints)
forceUpdate = bool(o.forceUpdate)

# Common parameters
isBump_win_dt = 125.
isBump_tstart = 0.
isBump_tend   = None
isBump_readme = 'Bump position estimation. Whole simulation'


# Create visitors
if common.bumpType in o.type:
    bumpVisitor = vis.bumps.BumpFittingVisitor(
        forceUpdate=forceUpdate,
        tstart='full',
        readme='Bump fitting. Whole simulation, starting at the start of theta stimulation.',
        bumpERoot='bump_e_full',
        bumpIRoot='bump_i_full')
    FRVisitor = vis.spikes.FiringRateVisitor(winLen=2.,     # ms
                                             winDt=.5,      # ms
                                             forceUpdate=forceUpdate)
    FRPlotter = vis.plotting.spikes.FiringRatePlotter(rootDir='pop_fr_plots')
    isBumpVisitor = vis.bumps.BumpPositionVisitor(
            tstart=isBump_tstart,
            tend=isBump_tend,
            win_dt=isBump_win_dt,
            readme=isBump_readme,
            forceUpdate=forceUpdate)

    sp.visit(bumpVisitor)
    sp.visit(isBumpVisitor)
    sp.visit(FRVisitor)
    #sp.visit(FRPlotter)

if common.gammaType in o.type:
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']

    statsVisitor_e = vis.spikes.SpikeStatsVisitor("spikeMon_e",
                                                  forceUpdate=forceUpdate)
    ACVisitor = vis.signals.AutoCorrelationVisitor(monName, stateList,
                                                   forceUpdate=forceUpdate)

    sp.visit(ACVisitor)
    sp.visit(statsVisitor_e)

if common.velocityType in o.type:
    speedEstimator = vis.bumps.SpeedEstimator(
            forceUpdate=forceUpdate,
            axis='vertical',
            win_dt=50.0)
    gainEstimator = vis.bumps.VelocityGainEstimator(
            o.bumpSpeedMax,
            forceUpdate=forceUpdate,
            maxFitRangeIdx=10)
    speedPlotter = vis.bumps.SpeedPlotter(plotFittedLine=True)

    sp.visit(speedEstimator, trialList='all-at-once')
    sp.visit(gainEstimator, trialList='all-at-once')
    sp.visit(speedPlotter, trialList='all-at-once')

if common.gridsType in o.type:

    po = vis.plotting.grids.GridPlotVisitor.PlotOptions()
    gridVisitor = vis.plotting.grids.GridPlotVisitor(o.output_dir,
                                                     spikeType='E',
                                                     plotOptions=po,
                                                     minGridnessT=300e3,
                                                     forceUpdate=o.forceUpdate)
    gridVisitor_i = vis.plotting.grids.IGridPlotVisitor(o.output_dir,
                                                        plotOptions=po,
                                                        minGridnessT=300e3,
                                                        forceUpdate=o.forceUpdate)
    isBumpVisitor = vis.bumps.BumpPositionVisitor(tstart=isBump_tstart,
                                                  tend=isBump_tend,
                                                  win_dt=isBump_win_dt,
                                                  readme=isBump_readme,
                                                  forceUpdate=forceUpdate,
                                                  bumpERoot='bump_e_isBump')
    #ISIVisitor = plotting_visitors.ISIPlotVisitor(o.output_dir,
    #        spikeType = spikeType,
    #        nRows = 5, nCols = 5, range=[0, 1000], bins=40,
    #        ISINWindows=20)
    FRVisitor = vis.spikes.FiringRateVisitor(winLen=2.,     # ms
                                             winDt=.5,      # ms
                                             forceUpdate=forceUpdate,
                                             sliding_analysis=False)

    sp.visit(gridVisitor)
    sp.visit(gridVisitor_i)
    #sp.visit(isBumpVisitor)
    #sp.visit(ISIVisitor)
    sp.visit(FRVisitor)

if common.gridsIPCType in o.type:
    # This is solely for the purpose of analyzing simulations where a
    # population of place cells is connected to I cells.

    po = vis.plotting.grids_ipc.GridPlotVisitor.PlotOptions()
    ipc_gridVisitor = vis.plotting.grids_ipc.GridPlotVisitor(
        o.output_dir,
        spikeType='E',
        plotOptions=po,
        minGridnessT=300e3,
        forceUpdate=o.forceUpdate)
    ipc_gridVisitor_i = vis.plotting.grids_ipc.IGridPlotVisitor(
        o.output_dir,
        plotOptions=po,
        minGridnessT=300e3,
        forceUpdate=o.forceUpdate)
    ipc_FRVisitor = vis.spikes.FiringRateVisitor(
        winLen=2.,     # ms
        winDt=.5,      # ms
        forceUpdate=forceUpdate,
        sliding_analysis=False)

    sp.visit(ipc_gridVisitor)
    sp.visit(ipc_gridVisitor_i)
    sp.visit(ipc_FRVisitor)

if common.posType in o.type:
    bumpPosVisitor = vis.bumps.BumpPositionVisitor(
            tstart=0,
            tend=None,
            win_dt=125.0,
            readme='Bump position estimation. Whole simulation.',
            forceUpdate=forceUpdate)
    sp.visit(bumpPosVisitor)


print('Total time: %.3f s' % (time.time() - startT))
