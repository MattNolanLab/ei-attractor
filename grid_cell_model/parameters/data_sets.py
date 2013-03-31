#
#   data_sets.py
#
#   Data set classes.
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


class DataSet(object):
    '''
    The data set object interface
    '''

    def __init__(self):
        raise NotImplementedError()

    @property
    def data(self):
        raise NotImplementedError()

    @property
    def parameters(self):
        raise NotImplementedError()



class DictDataSet(DataSet):

    def __init__(self, dataDict):
        self._d = dataDict


    @property
    def data(self):
        return self._d
