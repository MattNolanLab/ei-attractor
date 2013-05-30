#!/usr/bin/env python
#
#   details_EI.py
#
# Detailed Vm, Isyn, FR, and annotations plot for the 2D (EI) parameter sweep.
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

from analysis.visitors    import DetailedPlotVisitor
from parameters           import JobTrialSpace2D
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


###############################################################################
optParser = OptionParser()
optParser.add_option('--row',   type="int")
optParser.add_option('--col',   type="int")
optParser.add_option('--shapeRows',   type="int")
optParser.add_option('--shapeCols',   type="int")
optParser.add_option("--output_dir",  type="string")
optParser.add_option("--job_num",     type="int") # unused

o, args = optParser.parse_args()

###############################################################################

outRootDir = "{0}/{1}".format(o.output_dir, "detailed_plots")
plotT = .5e3 # ms
shape = (o.shapeRows, o.shapeCols)
dataPoints = [(o.row, o.col)]
trialNums = None

sp = JobTrialSpace2D(shape, o.output_dir, dataPoints=dataPoints)
detailVisitor = DetailedPlotVisitor(outRootDir, plotT)

sp.visit(detailVisitor, trialNums)
    
