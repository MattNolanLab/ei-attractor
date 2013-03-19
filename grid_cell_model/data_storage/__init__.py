#
#   __init__.py
#
#   Package initialization.
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


