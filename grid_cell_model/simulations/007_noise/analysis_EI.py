#!/usr/bin/env python
#
#   analysis_peaks.py
#
#   Theta/gamma analysis using a custom "peak" method.
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

from parameters import JobTrialSpace2D


###############################################################################

#rootDir = 'output_local/2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample'
rootDir = 'output/2013-03-30T19-34-44_EI_param_sweep_0pA_full/'
sp = JobTrialSpace2D((20, 20), rootDir)

