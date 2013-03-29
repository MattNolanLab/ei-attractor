#! /usr/bin/env python
#
#   test_submitting.py
#
#   Test the submitting package.
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
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator


ENV = 'workstation'
simRootDir =  'output_local'
simLabel   =  'test spaces'


parameters = {}


startJobNum = 10
numRepeat = 1


# Submitting
ac = ArgumentCreator(parameters)
appName    = 'dummy_script.py'
rtLimit    = '00:05:00'
blocking   = True
timePrefix = False
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix);
ac.setOption('output_dir', submitter.outputDir())


iterparams = {
        'output_data' : np.arange(10)
}
ac.insertDict(iterparams, mult=False)


submitter.submitAll(startJobNum, numRepeat, dry_run=False)
