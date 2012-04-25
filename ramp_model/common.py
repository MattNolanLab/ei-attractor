#
#   common.py
#
#     Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
'''This small module exports basic, shared classes and data'''

__all__ = ['NotImplementedException', 'DimensionException', 'ArgumentCreator',
        'ProgramSubmitter', 'GenericSubmitter', 'QsubSubmitter']

import os
from log import *



################################################################################
#                               Exceptions
################################################################################

class NotImplementedException(Exception):
    def __init__(self, method, msg=None):
        self.method = method
        self.msg = msg

    def __str__(self):
        retstr = "Method {0} has not been implemented!".format(self.method)
        if msg is not None:
            retstr += "Additional message: " + self.msg
        return retstr

class DimensionException(Exception):
    def __init__(self, method, msg=None):
        self.method = method
        self.msg = msg

    def __str__(self):
        retstr = "Wrong dimensions in {0}".format(self.method)
        if msg is not None:
            retstr += "Additional message: " + self.msg
        return retstr


################################################################################
#                      Argument processing and setup -
#             to submit jobs with multiple parameter sweeps/runs
################################################################################

class ArgumentCreator(object):
    '''
    Receives a dictionary with arguments and extracts a list of argument
    dictionaries, for submitting batch jobs. Either on a cluster or a general
    multiprocessor machine.

    The constructor receives a dictionary with default, non-batch changing
    optioins. After that, one can insert one of the following (see doc for each
    method)
    '''
    
    def __init__(self, defaultOpts):
        '''
        defaultOpts should be a dictionary (not a list)
        '''
        self._do = defaultOpts
        self._resList = [defaultOpts]

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
        with the value from 'd'. In other words, this allow to create
        multidimensional parameters batch jobs. If mult==false, then the size of
        list in d must be the same size as the total number of items in the
        resulting parameter list.
        '''
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


    def getOptionDict(self, i):
        '''
        Get option dictionary from the batch list, with index i
        '''
        return self._resList[i]

    def listSize(self):
        '''Return the size of the list of option dictionaries'''
        return len(self._resList)

    def getArgString(self, i):
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
        return parStr



class ProgramSubmitter(object):
    '''
    Submits arguments to the cluster or just a general machine. One should
    inherit from this class and implement/override the necessary methods.
    '''

    def __init__(self, argCreator):
        self._ac = argCreator

    def RunProgram(self, args):
        '''Should be overridden; args: argument string'''
        raise NotImplementedException("ProgramSubmitter.getProgramCommand")


    def submitAll(self):
        for it in range(self._ac.listSize()):
            print "Submitting simulation " + str(it)
            self.RunProgram(self._ac.getArgString(it))

    def submitOne(self, it):
        print "Submitting simulation " + str(it)
        self.RunProgram(self._ac.getArgString(it))


class GenericSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a generic multiprocessor machine, either by blocking them, or
    concurrently.
    '''
    def __init__(self, argCreator, appName, blocking=True):
        ProgramSubmitter.__init__(self, argCreator)
        self._appName = appName
        self._blocking = blocking

    def RunProgram(self, args):
        if self._blocking:
            postfix = ''
        else:
            postfix = '&'

        cmdStr = self._appName + ' ' + args + postfix
        log_info('root.GenericSubmitter', 'Submitting command: \n' + cmdStr)
        os.system(cmdStr)


class QsubSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a machine that supports qsub command (cluster)
    '''
    def __init__(self, argCreator, scriptName, qsub_params):
        ProgramSubmitter.__init__(self, argCreator)
        self._scriptName = scriptName
        self._qsub_params = qsub_params

    def RunProgram(self, args):
        cmd = 'qsub ' + qsub_params + ' '
        os.system(cmd + args)

