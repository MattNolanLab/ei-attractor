#!/usr/bin/env python
#
#   list_conversion.py
#
#   Conversion of list HDF5 data structures from the 0-prefix format to the
#   no-prefix format.
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
import h5py

def isGroup(ds):
    if (isinstance(ds, h5py.Group)):
        return True
    else:
        return False

def isList(grp):
    if (isinstance(grp, h5py.Group) and ('type' in grp.attrs.keys()) and
            grp.attrs['type'] == 'list'):
        return True
    else:
        return False


def fixList(grp):
    '''
    Convert all the indexes in group `grp`. From the 0-prefix format to a
    no-prefix format.
    '''
    print "  Fixing list (group) at {0}".format(grp.name)
    for oldIdx in grp.keys():
        newIdx = str(int(oldIdx))
        if (newIdx != oldIdx):
            grp[newIdx] = grp[oldIdx]
            del grp[oldIdx]

def fixGroup(grp):
    if (isGroup(grp)):
        if (isList(grp)):
            fixList(grp)
        # Even if the list was fixed, fix all of its children
        for key in grp.keys():
            fixGroup(grp[key])
    else:
        return


###############################################################################
if (len(sys.argv) != 2):
    exit(1)

fileList = list(open(sys.argv[1], 'r'))

it = 0
for fileName in fileList:
    print it
    if (fileName[-1] == '\n'):
        fileName = fileName[0:-1]
    print "Processing file: {0}".format(fileName)
    f = h5py.File(fileName, 'r+')
    fixGroup(f)
    f.close()
    it += 1

