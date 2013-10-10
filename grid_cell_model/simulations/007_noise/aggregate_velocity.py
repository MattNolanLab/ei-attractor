#!/usr/bin/env python
#
#   aggregate_velocity.py
#
#   Aggregate velocity analysis outputs.
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
import numpy as np
from parameters  import JobTrialSpace2D
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


noise_sigmas = [0, 150, 300]
dirs         = ('{0}pA', (31, 31))
shape        = dirs[1]
varListBase  = ['analysis']
loadData     = False

################################################################################
for noise_sigma in noise_sigmas:
    dir = dirs[0].format(int(noise_sigma))
    rootDir = "output/even_spacing/velocity/{0}".format(dir)

    sp = JobTrialSpace2D(shape, rootDir)
    sp.aggregateData(varListBase + ['lineFitErr'], funReduce=np.sum,
            trialNumList='all-at-once', saveData=True, output_dtype='array',
            loadData=loadData)
    sp.aggregateData(varListBase + ['lineFitSlope'], funReduce=None,
            trialNumList='all-at-once', saveData=True, output_dtype='array',
            loadData=loadData)
    sp.aggregateData(varListBase + ['bumpVelAll'], funReduce=None,
            trialNumList='all-at-once', saveData=True, output_dtype='list',
            loadData=loadData)
