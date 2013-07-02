#! /usr/bin/env python

import h5py

print(h5py.version.version)
print(h5py.version.hdf5_version)
print(h5py.version.api_version)

nm = 'test_h5py.h5'

f1 = h5py.File(nm, 'w')
f1.close()

f2 = h5py.File(nm, 'r')

f3 = h5py.File(nm, 'a')
f3.close()

f2.close()
