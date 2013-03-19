#
#   hdf5.py
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

    We use h5py here, because of its very simple design.
    '''

    def __init__(self, fileObj, grp):
        self._group = grp
        self._file  = fileObj


    @classmethod
    def factory(cls, filePath):
        f = h5py.File(filePath, 'a')
        return cls(f, f['/'])

       
    def __setitem__(self, key, value):
        if (key in self._group):
            del self._group[key] # TODO: this is costly operation
        self._createDataMember(key, value, self._group)



    def __getitem__(self, key):
        val = self._group[key]
        if (isinstance(val, h5py.Group)):
            return HDF5DataStorage(self._file, val)
        else:
            return val.value


    def __delitem__(self, key):
        del self._group[key]


    def __len__(self):
        return len(self._group)

    def __iter__(self):
        return iter(self._group.keys())


    def close(self):
        self._file.close()

    def flush(self):
        self._file.flush()


    def _createDataMember(self, name, value, grp):
        '''
        Assumes group has been created for us and we just fill it
        '''
        if (isinstance(value, dict)):
            newGrp = grp.create_group(name)
            for k, v in value.iteritems():
                self._createDataMember(k, v, newGrp)
        else:
            try:
                grp.create_dataset(name=name, data= value, compression="gzip")
            except TypeError:
                grp.create_dataset(name=name, data= value)


if __name__ == '__main__':
    import numpy as np

    ds = HDF5DataStorage.factory('test.h5')

    ds['ahoj'] = 124
    ds['compound'] = {'comp1': np.arange(100), 'comp2': 1234, 'str':
            'StringStringString',
            'nestedDict' : {'n1': 1, 'n2' : 2}}
    ds['big'] = np.random.rand(1000000)
    ds['2d'] = np.random.rand(100, 100)

    print ds['ahoj']

    ds.close()

    ds = HDF5DataStorage.factory('test.h5')
    print ds['compound']
    print ds
    print ds['compound']['comp1']
    print ds['compound']['comp2']
    print ds['compound']['str']
    print ds.keys()

    ds['compound']['comp1'][0] = 10
    print ds['compound']['comp1']

    ds.close()



