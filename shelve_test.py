#
#   shelve_test.py
#
#   Test the shelve module.
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
import shelve
import pickle
import numpy as np

filename = "data.shelve"
d = shelve.open(filename, protocol=pickle.HIGHEST_PROTOCOL)
d.clear()

d["ahoj"] = np.arange(1000000)

d.sync()
d.close()


del d

d = shelve.open(filename)
print d.keys()

arr = d['ahoj']
print arr
print arr.dtype

d.close()
