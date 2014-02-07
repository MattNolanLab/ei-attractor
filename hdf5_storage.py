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
from collections import MutableMapping, MutableSequence
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
    current implementation, the semantics of the **[]** operator (__getitem__)
    for list is that of a **shallow copy** of the object from the file.
    Therefore the following construct does not work::

        >>> d['listKey'][0] = 10
    This will change only the copy of the object stored in the file.

    Also, the following might have unexpected results, if the list stored in
    the HDF5 file is a compound list::
        
        1. >>> tmp = d['listKey']  # tmp now contains a shallow copy of the list
        3. >>> d['anotherKey'] = tmp

    The reason for this is that tmp contains only the first level shallow copy
    of the listKey list, for the sake of efficiency. However, compound objects
    like dictionaries will be returned as references, which will not be
    processed correctly. This implies that replacing a list somewhere in the
    data hierarchy with a shallow copy of the list (or another compound
    structure) will fail miserably.


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
    @staticmethod
    def factory(filePath, mode):
        '''
        Create the HDF5MapStorage object from a file.
        '''
        f = h5py.File(filePath, mode)
        return HDF5MapStorage(f, f['/'])


    def __init__(self, fileObj, grp):
        '''
        Initialize the object.

        Do not use this directly. Use the factory classmethod.
        '''
        self._group = grp
        self._file  = fileObj


    def _getitem(self, dataSet):
        '''
        Return the correct implementation of the data set `dataSet`.
        Three options are currently possible:
            1. Dict-like data set (in fact a hdf5 group, witth type attribute
               set to 'dict')
            2. List-like data set (a group again, with type attribute set to
               'list')
            3. An atomic data (everything else)

        1. and 2. will return a reference, 3. will return an actual copy of the
           object.
        '''
        val = dataSet
        if (isinstance(val, h5py.Group)):
            if (val.attrs['type'] == 'dict'):
                return HDF5MapStorage(self._file, val)
            elif (val.attrs['type'] == 'list'):
                return HDF5ListStorage(self._file, val)
            else:
                raise Exception("Unknown type attribute encountered while " +
                        "parsing the get request. Please check whether your " +
                        "HDF5 file is in correct format")
        else:
            try:
                return val.value
            except ValueError:
                if (len(val) == 0):
                    return np.array([], dtype=val.dtype)
                else:
                    raise


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
            if (isinstance(value, MutableMapping)):
                newGrp = grp.create_group(name)
                newGrp.attrs['type'] = 'dict'
                for k, v in value.iteritems():
                    self._createDataMember(k, v, newGrp)
            elif (isinstance(value, MutableSequence)):
                newGrp = grp.create_group(name)
                newGrp.attrs['type'] = 'list'
                it = 0
                for v in value:
                    self._createDataMember(str(it), v, newGrp)
                    it += 1
            else:
                try:
                    grp.create_dataset(name=name, data=value, compression="gzip")
                except TypeError:
                    grp.create_dataset(name=name, data=value)
                except ValueError:
                    mxShape = tuple([None]*value.ndim)
                    chunks = tuple([1]*value.ndim)
                    grp.create_dataset(name=name, shape=value.shape,
                            data=value, maxshape=mxShape, chunks=chunks,
                            compression="gzip")
        except TypeError:
            print "Could not create a data member %s" % name
            raise


    def getItemChained(self, keyTuple):
        '''
        Return an item at the and of a chain of keys, defined in ``keyTuple``
        '''
        l = len(keyTuple)
        if (l == 0):
            raise ValueError('Cannot chain index with an empty list of keys.')
        elif (l == 1):
            return self[keyTuple[0]]
        else:
            return self[keyTuple[0]].getItemChained(keyTuple[1:])


    def close(self):
        '''
        Close the file associated with this object.

        This will inevitably invalidate all the HDF5MapStorage objects derived
        from the same open file and subsequent access to any data members that
        are not copies will throw an exception.
        '''
        self._file.close()


    def flush(self):
        '''
        Flush all the data to the file stream.
        '''
        self._file.flush()



class HDF5MapStorage(HDF5DataStorage, MutableMapping):
    '''
    Dictionary-like HDF5 DataStorage implementation
    '''

    def __init__(self, fileObj, grp):
        HDF5DataStorage.__init__(self, fileObj, grp)
       
    def __setitem__(self, key, value):
        if (key in self._group):
            del self._group[key] # TODO: this is costly operation
        self._createDataMember(key, value, self._group)

    def __getitem__(self, key):
        #print("__getitem__({0}), group: {1}".format(key, self._group.name))
        return self._getitem(self._group[key])

    def __delitem__(self, key):
        del self._group[key]

    def __len__(self):
        return len(self._group)

    def __iter__(self):
        return iter(self._group.keys())


    def setItemChained(self, keyTuple, value, overwriteLast=True):
        '''
        Set ``value`` into ``keyTuple[-1]``. ``keyTuple`` must contain only
        strings, specifying dictionary keys. The semantic of this method is the
        following:

         * ``value`` can be of any supported types. It will simply be assigned
           as the last item in ``keyTuple``. If it already exists, it will be
           overwritten if ``overwriteLast`` is ``True``.
         * If ``keyTuple[0:-1]`` don't exist, create all of them as
           dictionaries. If they exist, do not overwrite.
        '''
        l = len(keyTuple)
        if (l == 0):
            raise ValueError("keyTuple must contain at least one item")
        elif (l == 1):
            if (not isinstance(keyTuple[0], str)):
                raise TypeError('All the keys in the keyTuple list must be strings.')
            if keyTuple[0] not in self.keys() or overwriteLast:
                self[keyTuple[0]] = value
        else:
            firstKey = keyTuple[0]
            if (not isinstance(firstKey, str)):
                raise TypeError('All the keys in the keyTuple list must be strings.')
            if firstKey in self.keys():
                self[firstKey].setItemChained(keyTuple[1:], value,
                        overwriteLast=overwriteLast)
            else:
                self[firstKey] = {}
                self[firstKey].setItemChained(keyTuple[1:], value,
                        overwriteLast=overwriteLast) 





class HDF5ListStorage(HDF5DataStorage, MutableSequence):
    '''
    List-like HDF5 DataStorage implementation.
    '''
    def __init__(self, fileObj, grp):
        HDF5DataStorage.__init__(self, fileObj, grp)
       
    def __setitem__(self, index, value):
        if (isinstance(index, slice)):
            raise TypeError('Slicing is not supported!')
        index = str(index)
        if (index in self._group):
            del self._group[index] # TODO: this is costly operation
            self._createDataMember(index, value, self._group)
        else:
            raise IndexError

    def __getitem__(self, index):
        if (isinstance(index, slice)):
            raise TypeError('Slicing is not supported!')
        if (index < 0):
            index = len(self._group) + index
        return self._getitem(self._group[str(index)])

    def __delitem__(self, index):
        raise NotImplementedError('Deleting items in an HDF5 stored list is' +\
                ' not as easy as you think!')

    def insert(self, index, value):
        raise NotImplementedError('Inserting items in an HDF5 stored list is' +\
                ' not as easy as you think!')

    def append(self, value):
        index = len(self)
        self._createDataMember(str(index), value, self._group)


    def setItemChained(self, keyTuple, value):
        '''
        This method does not make sence in HDF5ListStorage. Its semantic is
        only to create a hiearchy of dictionaries, e.g. as in
        :meth:`HDF5MapStorage.setItemChained`.
        '''
        raise RuntimeError("setItemChained() cannot be used here.")


    def __len__(self):
        return len(self._group)

    def __iter__(self):
        idx = 0
        stop = len(self)
        while (idx < stop):
            yield self[idx]
            idx += 1

    def __repr__(self):
        ret = '['
        for idx in xrange(len(self) - 1):
            ret += str(self[idx]) + ', '
        ret += str(self[-1]) + ']'
        return ret

    def __eq__(self, other):
        if (len(self) != len(other)):
            return False
        for idx in xrange(len(self)):
            if (self[idx] != other[idx]):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

