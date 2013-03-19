#! /usr/bin/env python2.6
#
#   data_storage.py
#
#   Test the data_storage package
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
from data_storage import DataStorage


ds = DataStorage.open('test.h5')

ds['ahoj'] = 124
ds['compound'] = {'comp1': np.arange(100), 'comp2': 1234, 'str':
        'StringStringString',
        'nestedDict' : {'n1': 1, 'n2' : 2}}
ds['big'] = np.arange(10000000)
ds['2d'] = np.random.rand(100, 100)

print "ds['ahoj']:\t", ds['ahoj']

ds.close()

ds = DataStorage.open('test.h5', mode='a')
print "ds['compound']\t", ds['compound']
print "ds\t", ds
print "ds['compound']['comp1']\t", ds['compound']['comp1']
print "ds['compound']['comp2']\t", ds['compound']['comp2']
print "ds['compound']['str']\t", ds['compound']['str']
print "ds.keys()", ds.keys()

ds['compound']['comp1'][0] = 10
print "ds['compound']['comp1']", ds['compound']['comp1']
print ds['compound']['comp1'][...]

ds.close()

