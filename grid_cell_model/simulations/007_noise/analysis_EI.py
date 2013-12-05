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
from optparse import OptionParser
import matplotlib
matplotlib.use('agg')

import analysis.visitors as vis
from parameters import JobTrialSpace2D

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

###############################################################################
optParser = OptionParser()
optParser.add_option('--row',          type="int")
optParser.add_option('--col',          type="int")
optParser.add_option('--shapeRows',    type="int")
optParser.add_option('--shapeCols',    type="int")
optParser.add_option('--forceUpdate',  type="int")
optParser.add_option("--output_dir",   type="string")
optParser.add_option("--job_num",      type="int") # unused
optParser.add_option("--type",         type="choice", choices=['gamma-bump', 'velocity', 'grids'])
optParser.add_option("--bumpSpeedMax", type="float")

o, args = optParser.parse_args()

###############################################################################

shape = (o.shapeRows, o.shapeCols)
dataPoints = [(o.row, o.col)]
trialNums = None

sp = JobTrialSpace2D(shape, o.output_dir, dataPoints=dataPoints)
forceUpdate = bool(o.forceUpdate)

if (o.type == "gamma-bump"):
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']
    iterList  = ['g_AMPA_total', 'g_GABA_total']
    ACVisitor = vis.AutoCorrelationVisitor(monName, stateList,
            forceUpdate=forceUpdate)
    bumpVisitor = vis.BumpFittingVisitor(forceUpdate=forceUpdate)
    FRVisitor = vis.FiringRateVisitor(forceUpdate=forceUpdate)
    CCVisitor = vis.CrossCorrelationVisitor(monName, stateList,
            forceUpdate=forceUpdate)

    #sp.visit(ACVisitor)
    sp.visit(bumpVisitor)
    #sp.visit(FRVisitor)
    #sp.visit(CCVisitor)
elif (o.type == "velocity"):
    VelVisitor = vis.BumpVelocityVisitor(o.bumpSpeedMax, forceUpdate=forceUpdate, printSlope=True)
    sp.visit(VelVisitor, trialList='all-at-once')
elif (o.type == 'grids'):
    spikeType = 'E'
    po = vis.GridPlotVisitor.PlotOptions()
    gridVisitor = vis.GridPlotVisitor(o.output_dir, spikeType=spikeType,
            plotOptions=po, minGridnessT=300e3)
    ISIVisitor = vis.ISIPlotVisitor(o.output_dir,
            spikeType = spikeType,
            nRows = 5, nCols = 5, range=[0, 1000], bins=40,
            ISINWindows=20)
    FRVisitor = vis.FiringRateVisitor(forceUpdate=forceUpdate)

    #sp.visit(gridVisitor)
    #sp.visit(ISIVisitor)
    sp.visit(FRVisitor)
else:
    raise ValueError("Unknown analysis type option: {0}".format(o.type))
