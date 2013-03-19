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
    '''
    pass


#class ShelveDataStorage(DataStorage):
#    pass
#
#
#class MatDataStorage(DataStorage):
#    '''
#    An implementation of DataStorage for the Matlab file format.
#    
#    Uses loadmat and savemat.
#    '''
#    pass




# Test this module
if __name__ == "__main__":
    pass
