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


ns_all  = [0, 150, 300]
ns_none = [-100]
varListBase  = ['analysis']
loadData     = True
dirs = \
    ('output/even_spacing/velocity_vertical/{0}pA', (31, 31), ns_all)
    #('output/detailed_noise/velocity/EI-3_1', (31, 9), ns_all)
    #('output/detailed_noise/velocity/EI-1_3', (31, 9), ns_all)
    #('output/even_spacing/velocity/{0}pA', (31, 31), ns_all)

################################################################################
shape        = dirs[1]
noise_sigmas = dirs[2]
for noise_sigma in noise_sigmas:
    rootDir = dirs[0].format(int(noise_sigma))

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
