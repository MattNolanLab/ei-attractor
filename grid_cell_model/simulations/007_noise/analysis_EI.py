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

from analysis.visitors    import AutoCorrelationVisitor, BumpFittingVisitor, \
        FiringRateVisitor
from parameters           import JobTrialSpace2D
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

###############################################################################
optParser = OptionParser()
optParser.add_option('--row',         type="int")
optParser.add_option('--col',         type="int")
optParser.add_option('--shapeRows',   type="int")
optParser.add_option('--shapeCols',   type="int")
optParser.add_option('--forceUpdate', type="int")
optParser.add_option("--output_dir",  type="string")
optParser.add_option("--job_num",     type="int") # unused

o, args = optParser.parse_args()

###############################################################################

shape = (o.shapeRows, o.shapeCols)
dataPoints = [(o.row, o.col)]
trialNums = None

sp = JobTrialSpace2D(shape, o.output_dir, dataPoints=dataPoints)

monName   = 'stateMonF_e'
stateList = ['I_clamp_GABA_A']
iterList  = ['g_AMPA_total', 'g_GABA_total']
forceUpdate = bool(o.forceUpdate)
ACVisitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=forceUpdate)
bumpVisitor = BumpFittingVisitor(forceUpdate=forceUpdate)
FRVisitor = FiringRateVisitor(forceUpdate=forceUpdate)

sp.visit(ACVisitor)
sp.visit(bumpVisitor)
sp.visit(FRVisitor)
