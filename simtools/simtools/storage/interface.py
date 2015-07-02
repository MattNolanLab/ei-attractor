'''Immplementation independent data storage interface classes (base
class(es)).
'''
from __future__ import absolute_import, print_function, division


class DataStorage(object):
    '''
    Class for saving and loading data structures transparently.

    Use this interface to access the data manipulation routines. This should be
    independent of the underlying data format.

    Use the open @classmethod (DataStorage.open(filePath)) to create the right
    implementation for your file format. Currently HDF5 files are supported
    only.

    Note that the [] operator returns a **copy** of a non-compound data
    structure, while if the data is compound (dict-like) it will return a
    reference only! This is done for the sake of efficiency, as loading big
    chunks of data from a dictionary can be quite inefficient. The copy
    semantics of the non-compound data access can be changed in the future.
    '''
    @staticmethod
    def open(filePath, mode='w'):
        '''
        Given the file Path, return the DataStorage object with the correct
        implementation, inferred from the filePath extension.

        filePath
            Full path to the filename to open
        mode='w'
            File mode. The default is to overwrite the old file. Use 'a' to
            append.
        '''
        from . import hdf5_storage
        return hdf5_storage.HDF5DataStorage.factory(filePath, mode=mode)
