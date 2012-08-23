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


import os
from log import *


# Get a preferred direction for a neuron
def getPreferredDirection(pos_x, pos_y):
# pos_x/y - position of neuron in 2d sheet
    pos4_x = pos_x % 2
    pos2_y = pos_y % 2
    if pos4_x == 0:
        if pos2_y == 0:
            return [-1, 0] # Left
        else:
            return [0, -1] # Down
    else:
        if pos2_y == 0:
            return [0, 1] # up
        else:
            return [1, 0] # Right



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

class SubmitException(Exception):
    pass

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
        self._printData = []

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

    def insertDict(self, d, mult=True, printout=False):
        '''
        Insert a dictionary with the parameter name and a list of values. If
        mult==True, then for each parameter value in 'd', copy the parameter set
        with the value from 'd'. In other words, this allow to create
        multidimensional parameters batch jobs. If mult==false, then the size of
        list in d must be the same size as the total number of items in the
        resulting parameter list.
        '''
        if printout:
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
            res += arg + ': {0:3.3f}'.format(self._resList[it][arg]) + '\t'
        return res



class ProgramSubmitter(object):
    '''
    Submits arguments to the cluster or just a general machine. One should
    inherit from this class and implement/override the necessary methods.
    '''

    def __init__(self, argCreator):
        self._ac = argCreator

    def RunProgram(self, args, job_num):
        '''Should be overridden; args: argument string'''
        raise NotImplementedException("ProgramSubmitter.getProgramCommand")


    def submitAll(self, startJobNum, repeat=1, dry_run=False):
        '''
        Submits all the generated jobs. Parameters:
          startJobNum   Start job number index
          repeat        Number of repeats for each parameters dictionary
        '''
        prt = []
        curr_job_num = startJobNum
        for it in range(self._ac.listSize()):
            for rep in range(repeat):
                print "Submitting simulation " + str(it)
                self.RunProgram(self._ac.getArgString(it, curr_job_num ), curr_job_num, dry_run)
                prt.append((curr_job_num, self._ac.getPrintArgString(it)))
                curr_job_num += 1

        for vals in prt:
            print vals[0], ': ',  vals[1]

    def submitOne(self, it, startJobNum, repeat=1):
        curr_job_num = startJobNum
        for rep in range(repeat):
            print "Submitting simulation " + str(it)
            self.RunProgram(self._ac.getArgString(it, curr_job_num), curr_job_num, dry_run)
            curr_job_num += 1


class GenericSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a generic multiprocessor machine, either by blocking them, or
    concurrently.
    '''
    def __init__(self, argCreator, appName, blocking=True):
        ProgramSubmitter.__init__(self, argCreator)
        self._appName = appName
        self._blocking = blocking

    def RunProgram(self, args, job_num, dry_run):
        if self._blocking:
            postfix = ''
        else:
            postfix = '&'

        cmdStr = self._appName + ' ' + args + postfix
        log_info('root.GenericSubmitter', 'Submitting command: \n' + cmdStr)
        if not dry_run:
            err = os.system(cmdStr)
            if err != 0:
                raise SubmitException()


class QsubSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a machine that supports qsub command (cluster)
    '''
    def __init__(self, argCreator, scriptName, qsub_params, qsub_output_dir):
        ProgramSubmitter.__init__(self, argCreator)
        self._scriptName        = scriptName
        self._qsub_params       = qsub_params
        self._qsub_output_dir   = qsub_output_dir

    def RunProgram(self, args, job_num, dry_run):
        cmdStr = 'qsub ' + self._qsub_params + ' ' + \
                "-o " + self._qsub_output_dir + ' ' + \
                "-N job{0:04} ".format(job_num) + \
                self._scriptName + ' ' + args
        log_info('root.QsubSubmitter', 'Submitting command: \n' + cmdStr)
        if not dry_run:
            err = os.system(cmdStr)
            if err != 0:
                raise SubmitException()

