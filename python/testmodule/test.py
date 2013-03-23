#
#   test.py
#
#   Simple test of the python-numpy_boost interface.
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
import _test


# int arrays
try:
    a1 = np.array([1, 2, 3, 4])
    a2 = np.array([1, 2, 3, 4])
    
    print "int:",  _test.sum(a1, a2)
except:
    print "Int arrays: Unexpected exception!"

# Double arrays
try:
    a1 = np.array([1., 2, 3, 4])
    a2 = np.array([1, 2, 3, 4])
    
    print "double:", _test.sum(a1, a2)
except:
    print "double arrays: Unexpected exception!"


try:
    a_short = np.array([1, 2, 3])
    _test.sum(a1, a_short)
    print "Different sizes: Failed! TypeError exception expected!"
except TypeError:
    print "Different sizes: OK: TypeError raised"


try:
    a1_tuple = (1, 2, 3)
    a2_tuple = (1, 2, 3)
    print "tuple:", _test.sum(a1_tuple, a2_tuple)
except:
    print "Tuple: Unexpected exception!"

