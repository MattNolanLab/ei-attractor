#
#   arguments.py
#
#   Argument manipulation.
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
'''
Argument processing and setup - to submit jobs with multiple parameter
sweeps/runs
'''

class ArgumentCreator(object):
    '''
    Receives a dictionary with arguments and extracts a list of argument
    dictionaries, for submitting batch jobs. Either on a cluster or a general
    multiprocessor machine.

    The constructor receives a dictionary with default, non-batch changing
    optioins. After that, one can insert one of the following (see doc for each
    method)
    '''
    
    def __init__(self, defaultOpts, printout=False):
        '''
        defaultOpts should be a dictionary (not a list)
        '''
        self._do = defaultOpts
        self._resList = [defaultOpts]
        self._printData = []
        self.printout = printout

        # remove job_num parameter from the list, it should be defined during
        # job submission
        try:
            del self._do["job_num"]
        except KeyError:
            pass

    def __str__(self):
        ret = ''
        for it in range(len(self._resList)):
            ret += 'Parameter set no. ' + str(it) + '\n'
            for key, val in self._resList[it].iteritems():
                ret += '\t' + key + ': ' + str(val) + '\n'

        return ret

    def insertDict(self, d, mult=True):
        '''
        Insert a dictionary with the parameter name and a list of values. If
        mult==True, then for each parameter value in 'd', copy the parameter set
        with the value from 'd'. In other words, this allows to create
        multidimensional parameters batch jobs. If mult==false, then the size of
        list in d must be the same size as the total number of items in the
        resulting parameter list.
        '''
        if (self.printout):
            for key, vals in d.iteritems():
                self._printData.append(key)

        if mult == True:
            # Create potentially multidim list (but flattened)
            for key, vals in d.iteritems():
                newList = []
                for val in vals:
                    for oldDict in self._resList:
                        newDict = dict(oldDict)
                        newDict[key] = val
                        newList.append(newDict)
                self._resList = newList
        else:
            # Check if list in d is the same size as _resList and insert
            if self.listSize() == 1:
                for it in range(len(d[d.keys()[0]]) - 1):
                    self._resList.append(dict(self._resList[0]))

            for key, vals in d.iteritems():
                if len(vals) != self.listSize():
                    raise DimensionException("ArgumentCreator.insertDict",
                            "Dictionary list must be the same size as the total number of items in the argument list")
                for it in range(len(vals)):
                    self._resList[it][key] = vals[it]

    def setOption(self, key, val):
        for opts in self._resList:
            opts[key] = val


    def getOptionDict(self, i):
        '''
        Get option dictionary from the batch list, with index i
        '''
        return self._resList[i]

    def listSize(self):
        '''Return the size of the list of option dictionaries'''
        return len(self._resList)

    def getArgString(self, i, job_num):
        '''
        Get the argument string which should be passed on to the submitted
        program (simulation).
        '''
        parStr = ''
        for key, val in self._resList[i].iteritems():
            if val is None:
                parStr += ' --' + str(key)
            else:
                parStr += ' --' + str(key) + ' ' + str(val)
        parStr += ' --job_num ' + str(job_num)
        return parStr

    def getPrintData(self):
        return self._printData

    def getPrintArgString(self, it):
        res = ''
        for arg in self._printData:
            res += arg + ' {0:3.3f}'.format(self._resList[it][arg]) + ' '
        return res



