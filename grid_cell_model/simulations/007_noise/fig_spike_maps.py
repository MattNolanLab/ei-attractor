#!/usr/bin/env python
#
#   fig_spike_maps.py
#
#   Spike, rate maps, etc. of grid cells.
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
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import rcParams

from parameters  import JobTrialSpace2D
from analysis.visitors import GridPlotVisitor
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

rcParams['font.size'] = 14
spikeType = 'E'

noise_sigma = 150
dirs = \
    ('grids_no_velocity_{0}pA',    (1, 100))

NTrials = 1
simDir = dirs[0].format(int(noise_sigma))
rootDir = "output/grids_no_velocity/{0}".format(simDir)
shape   = dirs[1]

sp = JobTrialSpace2D(shape, rootDir)#, dataPoints=[(0, 0)])
visitor = GridPlotVisitor(rootDir, spikeType=spikeType)
sp.visit(visitor)
