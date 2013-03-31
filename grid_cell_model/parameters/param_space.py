#
#   param_space.py
#
#   Parameter spaces classes.
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
from collections    import Sequence

from data_storage   import DataStorage
from data_sets      import DictDataSet

__all__ = []



class DataSpace(Sequence):
    '''
    An interface to the space of data with different parameters.
    '''

    def __init__(self, vals, pList=None, key=None):
        self._key  = key
        self._vals = vals
        self._pList = pList


    def __getitem__(self, key):
        return self._vals[key]


    def applyVisitor(self, visitor):
        for val in self._vals:
            val.applyVisitor(visitor)

    def __len__(self):
        return len(self._vals)



class TrialSet(DataSpace):
    '''
    A 1D DataSpace that contains a list of DataSet objects.
    '''

    def __init__(self, fileName):
        try:
            self._ds = DataStorage.open(fileName, 'r')
            DataSpace.__init__(self, self._ds['trials'], key='trials')
            print("Opened " + fileName)
        except IOError:
            DataSpace.__init__(self, [], key='trials')

    def __getitem__(self, key):
        return DictDataSet(self._vals[key])




class JobTrialSpace2D(DataSpace):
    def __init__(self, shape, rootDir, fileFormat="job{0:05}_output.h5"):
        self._shape = shape
        self._rootDir = rootDir
        self._fileFormat = fileFormat
        rows = shape[0]
        cols = shape[1]
        rowData = []
        it = 0
        for row in xrange(rows):
            colData = []
            for col in xrange(cols):
                fileName = rootDir + '/' + fileFormat.format(it)
                colData.append(TrialSet(fileName))
                it += 1
            rowData.append(DataSpace(colData))
        DataSpace.__init__(self, rowData)
            
                
    def __len__(self):
        return self.shape[0] 


