#
#   interface.py
#
#   Loading and saving data in the simulations/analysis.
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

from collections import MutableMapping


class DataStorage(MutableMapping):
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

    @classmethod
    def open(cls, filePath, mode='w'):
        '''
        Given the file Path, return the DataStorage object with the correct
        implementation, inferred from the filePath extension.
        
        filePath
            Full path to the filename to open
        mode='w'
            File mode. The default is to overwrite the old file. Use 'a' to
            append.
        '''
        import hdf5_storage
        return hdf5_storage.HDF5DataStorage.factory(filePath, mode=mode)




###############################################################################
#                                  Tests
###############################################################################
if __name__ == "__main__":
    pass
