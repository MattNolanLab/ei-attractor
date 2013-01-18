#
#   fig_param_sweep.py
#
#   Parameter sweep: data analysis and figures
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
from scipy      import linspace
from scipy.io   import loadmat

from tools                   import *
from plotting                import *

import numpy as np

output_dir = 'output'
fileNamePrefix = ''
job_num = 0
trial_it = 0


input_fname = "{0}/{1}job{2:05}_trial{3:04}_output.mat".format(output_dir,
        fileNamePrefix, job_num, trial_it)
inData = loadmat(input_fname)


