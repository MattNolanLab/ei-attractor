#! /usr/bin/env python
#
#   submit_bump_fitting.py
#
#   Submit job(s) to the cluster/workstation
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
from default_params         import defaultParameters as p
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator


# Submitting
ENV         = 'workstation'
simLabel    = 'tmp_bump_fitting'
simRootDir  = 'output_local'
appName     = 'simulation_stationary.py'
rtLimit     = '00:05:00'
blocking    = True
timePrefix  = False
numRepeat = 1



#p['time']              = 1199.9e3  # ms
p['time']              = 10e3  # ms
p['ntrials']           = 1
p['nthreads']          = 8

p['Ne']                = 34
p['Ni']                = 34

p['bumpCurrentSlope']  = 1.175     # pA/(cm/s), !! this will depend on prefDirC !!
p['gridSep']           = 70.0      # cm, grid field inter-peak distance
p['N_place_cells']     = 30*30

# Gamma analysis params
p['gammaNSample']      = 0.10      # fraction

p['noise_sigma']       = 150.0       # pA
p['delay']             = 0.1




###############################################################################
ac = ArgumentCreator(p)
iterp = {
        'noise_sigma' : [0.0, 150.0, 200, 250, 300]    # pA
        #'noise_sigma' : [0.0]
}
ac.insertDict(iterp, mult=False)


###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix);
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
