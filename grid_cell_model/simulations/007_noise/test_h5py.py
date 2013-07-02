#! /usr/bin/env python

import h5py

print(h5py.version.version)
print(h5py.version.hdf5_version)
print(h5py.version.api_version)


nm = 'test_h5py.h5'
try:
    f = h5py.File(nm, 'r')
    var = f['non-existent']
    f.close()
except:
    print("Exception happened")
    #f.close()


f1 = h5py.File(nm, 'a')
f1.close()
