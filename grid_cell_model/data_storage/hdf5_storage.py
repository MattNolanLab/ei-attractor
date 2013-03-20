#
#   hdf5_storage.py
#
#   Data storage using the HDF5 library (and h5py currently)
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
import h5py
from interface import DataStorage


class HDF5DataStorage(DataStorage):
    '''
    An implementation of DataStorage for the HDF5 data format.

    Provides a dictionary interface to the HDF5 backend, which is currently
    implemented using h5py, as it maps very nicely to the dictionary-like
    structure.

    These types of objects can be stored currently:
        * primitive numeric data types
        * array-like objects (ideally as numpy arrays)
        * lists (see note on `List performance`_
        * strings
        * Objects of class dict and classes derived from dict, containing all
            of the above types, including itself

    All objects of non-dict type are stored in datasets, while for each dict a
    group is created (with a corresponding name).

    Note that, while the interface is the same as that of a dictionary, in the
    current implementation, the semantics of a **[]** operator (__getitem__) is
    that of a **copy** of the object from the file. Therefore the following
    construct does not work::

        >>> d[key'][0] = 10

    This will change only the copy of the object stored in the file.

    List performance
    ----------------
    Lists in Python are heterogeneous data structure and therefore the storage
    in HDF5 is not trivial. Currently a list is stored as a group with an
    attribute 'type' == 'list'. Each item in a list can is stored in a
    corresponding format (i.e. dictionaries and lists as groups, other types
    natively). This, indeed, has some performance limitations: creation and
    reading of large lists of basic data types will be very slow (even simple
    lists of 5000 elements will take a few seconds to store and load).

    It is therefore recommended to convert the list to some other data type
    before storing it through this interface.

    Cycles
    ------
    Note that in the current version, cycles are not detected in compound
    objects(dict, list). The user must handle these situations beforhand.
    '''

    def __init__(self, fileObj, grp):
        '''
        Initialize the object.

        Do not use this directly. Use the factory classmethod.
        '''
        self._group = grp
        self._file  = fileObj


    @classmethod
    def factory(cls, filePath, mode):
        '''
        Create the HDF5DataStorage object from a file.
        '''
        f = h5py.File(filePath, mode)
        return cls(f, f['/'])

       
    def __setitem__(self, key, value):
        if (key in self._group):
            del self._group[key] # TODO: this is costly operation
        self._createDataMember(key, value, self._group)

    def __getitem__(self, key):
        val = self._group[key]
        if (isinstance(val, h5py.Group)):
            if (val.attrs['type'] == 'dict'):
                return HDF5DataStorage(self._file, val)
            elif (val.attrs['type'] == 'list'):
                res = []
                for k in val.keys():
                    res.append(HDF5DataStorage(self._file, val)[k])
                return res
            else:
                raise Exception("Unknown type attribute encountered while " +
                        "parsing the get request. Please check whether your " +
                        "HDF5 file is in correct format")
        else:
            return val.value

    def __delitem__(self, key):
        del self._group[key]

    def __len__(self):
        return len(self._group)

    def __iter__(self):
        return iter(self._group.keys())


    def close(self):
        '''
        Close the file associated with this object.

        This will inevitably invalidate all the HDF5DataStorage objects derived
        from the same open file and subsequent access to any data members that
        are not copies will throw an exception.
        '''
        self._file.close()

    def flush(self):
        '''
        Flush all the data to the file stream.
        '''
        self._file.flush()


    def _createDataMember(self, name, value, grp):
        '''
        Create data member with name, value, in group grp.

        This assumes that 'name' can be safely created in grp.

        We distinguish three data types:
            * non-compound: atomic types, arrays, strings
            * compound: dict objects and classes inherited from dict. The
                compound types are stored recursively and for each a group is
                created instead of dataset.
            * lists: list objects and classes inherited from list. A list is
              stored recursively as a group. See a note in `List performance`_
              for the performance limitations of storing lists.
        '''
        try:
            if (isinstance(value, dict)):
                newGrp = grp.create_group(name)
                newGrp.attrs['type'] = 'dict'
                for k, v in value.iteritems():
                    self._createDataMember(k, v, newGrp)
            elif (isinstance(value, list)):
                newGrp = grp.create_group(name)
                newGrp.attrs['type'] = 'list'
                ln = len(value)
                if (ln != 0):
                    nDigits = int(np.ceil(np.log10(ln)))
                    digit = "{0:0" + str(nDigits) + "}"
                it = 0
                for v in value:
                    self._createDataMember(digit.format(it), v, newGrp)
                    it += 1
            else:
                try:
                    grp.create_dataset(name=name, data= value, compression="gzip")
                except TypeError:
                    grp.create_dataset(name=name, data= value)
        except TypeError:
            print "Could not create a data member %s" % name
            raise


