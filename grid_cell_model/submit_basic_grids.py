#!/usr/bin/env python
#
#   submit_basic_grids.py
#
#   Submit job(s) to the cluster/workstation:
#     Grid fields with theta input and all the inhibition (for gamma) and place
#     input.
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
from common         import GenericSubmitter, ArgumentCreator

import logging as lg
lg.basicConfig(level=lg.DEBUG)


parameters = defaultParameters

#parameters['time']              = 1199.9e3  # ms
parameters['time']              = 10e3      # ms
parameters['ndumps']            = 1

parameters['placeT']            = 10e3      # ms

parameters['stateMonDur']       = 10e3

parameters['bumpCurrentSlope']  = 1.175     # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 60        # cm, grid field inter-peak distance

startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python simulation_basic_grids.py'
blocking            = True

ac = ArgumentCreator(parameters)
submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
