#! /usr/bin/env python

import h5py

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
