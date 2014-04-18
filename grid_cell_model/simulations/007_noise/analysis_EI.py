#!/usr/bin/env python
#
'''
Perform analysis on whole 2D data sets.
'''
import time
import matplotlib
matplotlib.use('agg')

import visitors as vis
import visitors.spikes
import visitors.bumps
import visitors.signals
import visitors.plotting
import visitors.plotting.spikes
import visitors.plotting.grids
from parameters import JobTrialSpace2D
from submitting import flagparse
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
parser.add_argument("--type",         type=str, choices=common.allowedTypes, required=True)
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
if o.type == common.bumpType:
    bumpVisitor = vis.bumps.BumpFittingVisitor(forceUpdate=forceUpdate,
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

    #sp.visit(bumpVisitor)
    sp.visit(isBumpVisitor)
    sp.visit(FRVisitor)
    sp.visit(FRPlotter)

elif o.type == common.gammaType:
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']

    statsVisitor_e = vis.spikes.SpikeStatsVisitor("spikeMon_e",
                                                  forceUpdate=forceUpdate)
    ACVisitor = vis.signals.AutoCorrelationVisitor(monName, stateList,
                                                   forceUpdate=forceUpdate)
    CCVisitor = vis.signals.CrossCorrelationVisitor(monName, stateList,
                                                    forceUpdate=forceUpdate)

    sp.visit(ACVisitor)
    sp.visit(CCVisitor)
    sp.visit(statsVisitor_e)

elif o.type == common.velocityType:
    speedEstimator = vis.bumps.SpeedEstimator(
            forceUpdate=forceUpdate,
            axis='vertical',
            win_dt=50.0)
    gainEstimator = vis.bumps.VelocityGainEstimator(
            o.bumpSpeedMax,
            forceUpdate=forceUpdate,
            maxFitRangeIdx=10)
    speedPlotter = vis.bumps.SpeedPlotter(plotFittedLine=False)

    sp.visit(speedEstimator, trialList='all-at-once')
    #sp.visit(gainEstimator, trialList='all-at-once')
    #sp.visit(speedPlotter, trialList='all-at-once')

elif o.type == common.gridsType:
    spikeType = 'E'

    po = vis.plotting.grids.GridPlotVisitor.PlotOptions()
    gridVisitor = vis.plotting.grids.GridPlotVisitor(
            o.output_dir,
            spikeType=spikeType,
            plotOptions=po,
            minGridnessT=300e3,
            forceUpdate=o.forceUpdate)
    isBumpVisitor = vis.bumps.BumpPositionVisitor(
            tstart=isBump_tstart,
            tend=isBump_tend,
            win_dt=isBump_win_dt,
            readme=isBump_readme,
            forceUpdate=forceUpdate,
            bumpERoot='bump_e_isBump')
    #ISIVisitor = plotting_visitors.ISIPlotVisitor(o.output_dir,
    #        spikeType = spikeType,
    #        nRows = 5, nCols = 5, range=[0, 1000], bins=40,
    #        ISINWindows=20)
    #FRVisitor = plotting_visitors.FiringRateVisitor(forceUpdate=forceUpdate)

    sp.visit(gridVisitor)
    sp.visit(isBumpVisitor)
    #sp.visit(ISIVisitor)
    #sp.visit(FRVisitor)
elif o.type == common.posType:
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
