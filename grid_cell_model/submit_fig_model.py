#!/usr/bin/env python
#
#   submit_fig_model.py
#
#   Submit job(s) to the cluster/workstation: basic model figures
#     Model description figures
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
from default_params import defaultParameters
from common         import *
import logging  as lg
import numpy    as np
lg.basicConfig(level=lg.DEBUG)


parameters = defaultParameters
parameters['time']              = 6e3   # ms
parameters['theta_start_mon_t'] = 1e3   # ms


startJobNum = 100
numRepeat = 1

# Workstation parameters
programName = 'python simulation_fig_model.py'
blocking    = True
ac = ArgumentCreator(parameters)
submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
