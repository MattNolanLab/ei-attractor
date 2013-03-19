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

        We distinguish two data types:
            * non-compound: atomic types, arrays, strings
            * compound: dict objects and classes inherited from dict. The
                compound types are stored recursively and for each a group is
                created instead of dataset.
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
                it = 0
                for v in value:
                    self._createDataMember(str(it), v, newGrp)
                    it += 1
            else:
                try:
                    grp.create_dataset(name=name, data= value, compression="gzip")
                except TypeError:
                    grp.create_dataset(name=name, data= value)
        except TypeError:
            print "Could not create a data member %s" % name
            raise


