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
from analysis.visitors    import DetailedPlotVisitor
from parameters           import JobTrialSpace2D
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


###############################################################################

dirs = {
    './output_local/2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample' : (2, 2),
    #'./output_local/2013-04-24T15-27-30_EI_param_sweep_0pA_big'   : (40, 40),
    #'./output_local/2013-04-24T21-37-47_EI_param_sweep_150pA_big' : (40, 40),
    #'output_local/2013-04-24T21-43-32_EI_param_sweep_300pA_big' : (40, 40)
}
NTrials = 5

for rootDir, shape in dirs.iteritems():
    outRootDir = "{0}/{1}".format(rootDir, "detailed_plots")
    plotT = .5e3 # ms
    #dataPoints = [(0, 0)]
    dataPoints = None
    #trialNums = [0]
    trialNums = None

    sp = JobTrialSpace2D(shape, rootDir, dataPoints=dataPoints)
    detailVisitor = DetailedPlotVisitor(outRootDir, plotT)
    
    sp.visit(detailVisitor, trialNums)
        
