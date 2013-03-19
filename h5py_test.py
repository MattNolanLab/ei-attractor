#
#   h5py_test.py
#
#   Test HDF5 access.
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


import h5py
import numpy as np

from collections import MutableMapping


fileName = "test_h5py.h5"
f = h5py.File(fileName, 'w')

grp = f.create_group('testGrp')
print grp.name
print grp.attrs.keys()

grp.attrs['attr1'] = 'ahoj'
print grp.attrs.keys()
print grp.attrs['attr1']

ds = f.create_dataset('dtst', data=np.array([10, 20.0, 30]))

print f['dtst']
print f['dtst'][:]
print f['dtst'].dtype
print f['dtst'].shape

print f.keys()
print len(f)

del f['dtst']
print len(f)

del f['testGrp']
print len(f)



grp = f.create_group("duplicate")
f['duplicate']['ahoj'] = 1234

print f['duplicate']['ahoj'].value
print f['duplicate']['ahoj']
print type(f['duplicate']['ahoj'].value)



a = {'a': 1, 'b': 2}
print a.__iter__()


isscalar(a)

