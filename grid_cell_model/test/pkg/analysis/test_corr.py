#! /usr/bin/env python
#
#   test_corr.py
#
#   Test the analysis.signal.corr function
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
import sys
import numpy as np

from analysis.signal import corr
from numpy           import correlate


def checkBlitzNumpyCorr(a1, a2, tolerance, printOut=False):
    c_blitz = corr(a1, a2, mode='twosided')
    c_np    = correlate(a1, a2, mode='full')[::-1]
    
    if (printOut):
        print c_blitz
        print c_np

    neq = c_blitz != c_np
    if (np.any(neq) and np.any(np.abs(c_blitz - c_np) > tolerance)):
        print c_blitz
        print c_np
        print c_blitz != c_np
        print c_blitz - c_np
        raise Exception("Arange: c_blitz != v_np!!")
###############################################################################

if (len(sys.argv) > 1):
    tolerance = float(sys.argv[1])
    print "Using tolerance= " + str(tolerance)
else:
    tolerance = .0

###############################################################################

N1 = 10
N2 = 6
a1 = np.arange(N1) * 1.0
a2 = np.arange(N2) * 1.0

checkBlitzNumpyCorr(a1, a2, tolerance)

###############################################################################
maxN = 1000
maxLoops = 10000
#maxLoops = np.infty

loop_it = 0
while (loop_it < maxLoops):
    N1 = np.random.randint(maxN) + 1
    N2 = np.random.randint(maxN) + 1
    if (N1 == 0 and N2 == 0):
        continue

    print N1, N2

    a1 = np.random.rand(N1)
    a2 = np.random.rand(N2)
    
    checkBlitzNumpyCorr(a1, a2, tolerance)

    loop_it += 1

