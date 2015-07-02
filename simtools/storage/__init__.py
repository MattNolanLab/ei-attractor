'''
====================
Data Storage package
====================

The data storage package contains interfaces and classes for storing simulation
output data independently of the file format.

Given that most of the time users generate large datasets (especially when
running large simulations), the target storage medium is a filesystem
directory (and in case of HDF it could be a single file).

Each simulation result is meant to be stored as a separate set of data (i.e. a
set of objects, each having its identifier). We will use the term *data set*
for this (not to be confused with HDF5 dataset concept).

A data set has a unique identifier, a string containing its name which is of
the form "name.extension". From this, the underlying data storage mechanism
will be inferred: for example, if a user wants to save the data as a HDF5 file,
they will use "name.h5" as an identifier. If the user wants to use the shelve
module, they will just use "name.shelve", etc.
'''
from __future__ import absolute_import, print_function, division

from .interface import DataStorage
