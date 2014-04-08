#!/usr/bin/env python
#
'''
Perform analysis on whole 2D data sets.
'''
import time
import matplotlib
matplotlib.use('agg')

import visitors as vis
import visitors.bumps as bumps
import visitors.spikes
import visitors.plotting.spikes
from parameters import JobTrialSpace2D
from submitting import flagparse

###############################################################################
parser = flagparse.FlagParser()
parser.add_argument('--row',          type=int, required=True)
parser.add_argument('--col',          type=int, required=True)
parser.add_argument('--shapeRows',    type=int, required=True)
parser.add_argument('--shapeCols',    type=int, required=True)
parser.add_argument('--forceUpdate',  type=int, required=True)
parser.add_argument("--output_dir",   type=str, required=True)
parser.add_argument("--job_num",      type=int) # unused
parser.add_argument("--type",         type=str,
        choices=['gamma-bump', 'velocity', 'grids', 'positional'], required=True)
parser.add_argument("--bumpSpeedMax", type=float)

o = parser.parse_args()

###############################################################################
startT = time.time()

shape = (o.shapeRows, o.shapeCols)
dataPoints = [(o.row, o.col)]
trialNums = None

sp = JobTrialSpace2D(shape, o.output_dir, dataPoints=dataPoints)
forceUpdate = bool(o.forceUpdate)

# Create visitors
if (o.type == "gamma-bump"):
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']
    iterList  = ['g_AMPA_total', 'g_GABA_total']
    #ACVisitor = vis.analysis_visitors.AutoCorrelationVisitor(monName, stateList,
    #        forceUpdate=forceUpdate)
    bumpVisitor = vis.bumps.BumpFittingVisitor(forceUpdate=forceUpdate,
            tstart='full',
            readme='Bump fitting. Whole simulation, starting at the start of theta stimulation.',
            bumpERoot='bump_e_full',
            bumpIRoot='bump_i_full')
    bumpPosVisitor = vis.bumps.BumpPositionVisitor(
            tstart=0,
            tend=None,
            win_dt=125.0,
            readme='Bump position estimation. Whole simulation',
            forceUpdate=forceUpdate)
    FRVisitor = vis.spikes.FiringRateVisitor(winLen=2.,     # ms
                                             winDt=.5,      # ms
                                             forceUpdate=forceUpdate)
    FRPlotter = vis.plotting.spikes.FiringRatePlotter(rootDir='pop_fr_plots')
    #CCVisitor = vis.CrossCorrelationVisitor(monName, stateList,
    #        forceUpdate=forceUpdate)
    #spikeVisitor_e = vis.SpikeStatsVisitor("spikeMon_e",
    #        forceUpdate=forceUpdate)

    #sp.visit(ACVisitor)
    #sp.visit(bumpVisitor)
    sp.visit(bumpPosVisitor)
    sp.visit(FRVisitor)
    sp.visit(FRPlotter)
    #sp.visit(CCVisitor)
    #sp.visit(spikeVisitor_e)
elif (o.type == "velocity"):
    speedEstimator = bumps.SpeedEstimator(
            forceUpdate=forceUpdate, axis='vertical', win_dt=50.0)
    gainEstimator = bumps.VelocityGainEstimator(o.bumpSpeedMax,
                                                forceUpdate=forceUpdate,
                                                maxFitRangeIdx=10)
    speedPlotter = bumps.SpeedPlotter(plotFittedLine=False)

    sp.visit(speedEstimator, trialList='all-at-once')
    #sp.visit(gainEstimator, trialList='all-at-once')
    #sp.visit(speedPlotter, trialList='all-at-once')
elif (o.type == 'grids'):
    spikeType = 'E'
    #po = plotting_visitors.GridPlotVisitor.PlotOptions()
    #gridVisitor = plotting_visitors.GridPlotVisitor(o.output_dir, spikeType=spikeType,
    #        plotOptions=po, minGridnessT=300e3, forceUpdate=o.forceUpdate)
    #ISIVisitor = plotting_visitors.ISIPlotVisitor(o.output_dir,
    #        spikeType = spikeType,
    #        nRows = 5, nCols = 5, range=[0, 1000], bins=40,
    #        ISINWindows=20)
    #FRVisitor = plotting_visitors.FiringRateVisitor(forceUpdate=forceUpdate)

    #sp.visit(gridVisitor)
    #sp.visit(ISIVisitor)
    #sp.visit(FRVisitor)
elif o.type == 'positional':
    bumpPosVisitor = vis.bumps.BumpPositionVisitor(
            tstart=0,
            tend=None,
            win_dt=125.0,
            readme='Bump position estimation. Whole simulation.',
            forceUpdate=forceUpdate)
    sp.visit(bumpPosVisitor)
else:
    raise ValueError("Unknown analysis type option: {0}".format(o.type))

print('Total time: %.3f s' % (time.time() - startT))
