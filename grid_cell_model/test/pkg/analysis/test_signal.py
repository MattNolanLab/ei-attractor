#
#   test_signal.py
#
#   Test the analysis.signal module
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
from analysis.signal import globalExtrema


###############################################################################
# globalExtrema
arr_1d = np.arange(100)
print globalExtrema(arr_1d, np.max)

Nrows = 11
Ncols = 100
arr_2d = np.ndarray((Nrows, Ncols))
arr_2d[:, :] = np.arange(100)
res_max = globalExtrema(arr_2d, np.max)
if (not np.all(res_max == Ncols - 1)):
    raise Exception("globalExtrema: 2d max failed!:")
else:
    print res_max

res_min = globalExtrema(arr_2d, np.min)
if (not np.all(res_min == 0)):
    raise Exception("globalExtrema: 2d min failed!:")
else:
    print res_min

arr_3d = np.ndarray((10, 10, 10))
try:
    globalExtrema(arr_3d)
    raise Exception("Test failed: globalExtrema 3d array!")
except TypeError:
    print "globalExtrema: 3d array: OK, TypError raised."


arr_empty = np.ndarray((0, ))
try:
    print globalExtrema(arr_empty, np.max)
    raise Exception("Test failed: globalExtrema empty array: expecting " +
            "ValueError exception")
except:
    print "globalExtrema: empty array: OK, ValueError raised."


