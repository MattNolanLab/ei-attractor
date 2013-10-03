#!/usr/bin/env python
#
#   submit_param_sweep_connections.py
#
#   Submit job(s) to the cluster/workstation: export connections to file
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
from default_params import defaultParameters as dp
from param_sweep    import submitParamSweep
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)


p = dp.copy()

# Submitting
ENV         = 'workstation'
simRootDir  = 'output_local/even_spacing'
simLabel    = 'connections'
appName     = 'simulation_connections.py'
rtLimit     = '00:02:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['nthreads'] = 1
p['ntrials']  = 1
p['connOnly'] = 1
#p['connNE']   = 10
#p['connNI']   = p['connNE']


# Range of E/I synaptic conductances
Nvals  = 31      # Number of values for each dimension
startG = 0.0     # nS
endG   = 6120.0  # nS

###############################################################################

submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
        appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run)
