#
#   visitors.py
#
#   Data analysis visitors. 
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


class Visitor(object):
    '''
    An abstract visitor class.

    Normally, the base class of the Visitor design pattern must contain all the
    methods implemented in the derived classes. Due to duck typing, one does
    not need to declare the specific implementation methods of the visitor.
    '''
    def __init__(self):
        raise NotImplementedError()



class DictDSVisitor(Visitor):
    '''
    Dictionary data set visitor.

    A visitor that takes a dictionary data set method of any kind and processes
    it. All the keys in the dictionary must be strings.
    '''
    def __init__(self):
        raise NotImplementedError()

    def visitDictDataSet(self, ds):
        '''
        Visit the dictionary data set, 'ds', and perform the specific operations
        (defined by the derived classes) on this data set.
        '''
        raise NotImplementedError()

    def getOption(self, data, optStr):
        '''Extract an option from a data dictionary'''
        return data['options'][optStr]

    def getNetParam(self, data, p):
        '''Extract a network parameter (p) from the data dictionary'''
        return data['net_attr'][p]


    def folderExists(self, d, nameList):
        '''
        Check if the folder at the end of the 'nameList' exists in dictionary
        'd'
        '''
        for name in nameList:
            if (name in d.keys()):
                d = d[name]
            else:
                return False
        return True

    def _checkAttrIsNone(self, attr, pName, data):
        '''
        Check if 'attr' is None. If yes, extract it from 'data' under the name
        'pName', otherwise return attr back.
        '''
        if (attr is None):
            return self.getOption(data, pName)
        else:
            return attr


