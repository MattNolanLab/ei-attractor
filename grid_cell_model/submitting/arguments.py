'''Argument processing and setup.

To submit jobs with multiple parameter sweeps/runs.
'''
from __future__ import absolute_import, print_function

from ..gc_exceptions import DimensionException

class ArgumentCreator(object):
    '''
    Receives a dictionary with arguments and extracts a list of argument
    dictionaries, for submitting batch jobs. Either on a cluster or a general
    multiprocessor machine.

    The constructor receives a dictionary with default, non-batch changing
    optioins. After that, one can insert one of the following (see doc for each
    method)
    '''
    def __init__(self, defaultOpts, printout=False, emitJobNum=True):
        '''
        defaultOpts should be a dictionary (not a list)
        '''
        self._do = defaultOpts
        self._resList = [defaultOpts]
        self._printData = []
        self.printout = printout
        self.emitJobNum = emitJobNum

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

    def optionList(self):
        '''Get the list of all options'''
        return self._resList

    def listSize(self):
        '''Return the size of the list of option dictionaries'''
        return len(self._resList)

    def _valToString(self, val):
        ret = ""
        if isinstance(val, list):
            for v in val:
                ret += str(v) + " "
        else:
            ret = str(val)
        return ret

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
                parStr += ' --' + str(key) + ' ' + self._valToString(val)
        if self.emitJobNum:
            parStr += ' --job_num ' + str(job_num)
        return parStr

    def getPrintData(self):
        return self._printData

    def getPrintArgString(self, it):
        res = ''
        for arg in self._printData:
            res += arg + ' {0:3.3f}'.format(self._resList[it][arg]) + ' '
        return res



