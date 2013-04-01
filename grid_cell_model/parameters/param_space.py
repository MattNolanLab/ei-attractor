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
from otherpkg.log   import log_warn, log_info

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


    def visit(self, visitor):
        for val in self._vals:
            val.visit(visitor)

    def __len__(self):
        return len(self._vals)



class TrialSet(DataSpace):
    '''
    A 1D DataSpace that contains a list of DataSet objects.
    '''

    def __init__(self, fileName, fileMode):
        self._fileName = fileName
        try:
            self._ds = DataStorage.open(fileName, fileMode)
            DataSpace.__init__(self, self._ds['trials'], key='trials')
            log_info("param_space", "Opened " + fileName)
        except IOError:
            self._ds = None
            msg =  "Could not open file {0}. Creating an empty DataSet instead."
            log_warn("param_space", msg.format(fileName))
            DataSpace.__init__(self, [], key='trials')

    def __del__(self):
        if (self._ds is not None):
            log_info("param_space", "Closing: " + self._fileName)
            self._ds.close()

    def __getitem__(self, key):
        return DictDataSet(self._vals[key])

    def visit(self, visitor):
        for val in self:
            val.visit(visitor)




class JobTrialSpace2D(DataSpace):
    def __init__(self, shape, rootDir, fileMode='r+',
            fileFormat="job{0:05}_output.h5", forceWMode=False):
        if (fileMode == 'w' and forceWMode == False):
            raise ValueError("'w' file open mode is not allowed. Use " +
                    "'forceWMode' to override.")
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
                colData.append(TrialSet(fileName, fileMode))
                it += 1
            rowData.append(DataSpace(colData))
        DataSpace.__init__(self, rowData)
            
                
    def __len__(self):
        return self.shape[0] 


    def getShape(self):
        return self._shape

