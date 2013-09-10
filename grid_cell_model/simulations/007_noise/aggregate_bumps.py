#!/usr/bin/env python
#
#   aggregate_bumps.py
#
#   Aggregate bump data into the reductions file.
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
from parameters  import JobTrialSpace2D
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


noise_sigmas = [150, 300]
dirs = \
    ('EI_param_sweep_{0}pA',    (40, 40))

NTrials = 5
trialNumList = xrange(NTrials)
shape   = dirs[1]
varListBase = ['analysis']

################################################################################
for noise_sigma in noise_sigmas:
    dir = dirs[0].format(int(noise_sigma))
    rootDir = "output/one_to_one/{0}".format(dir)

    sp = JobTrialSpace2D(shape, rootDir)
    sp.aggregateData(varListBase + ['bump_e', 'sigma'], trialNumList,
            funReduce=None, saveData=True, output_dtype='array')
    sp.aggregateData(varListBase + ['bump_e', 'err2'], trialNumList,
            funReduce=None, saveData=True, output_dtype='array')
    sp.aggregateData(varListBase + ['bump_e', 'bump_e_rateMap'], trialNumList,
            funReduce=None, saveData=True, output_dtype='list')
