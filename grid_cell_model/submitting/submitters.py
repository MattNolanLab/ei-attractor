from __future__ import absolute_import

import os, subprocess
import time
import errno
from datetime       import datetime

from simtools.storage import DataStorage

from ..gc_exceptions  import SubmitError
from ..otherpkg.log   import log_warn, log_info, getClassLogger


progSLogger = getClassLogger("ProgramSubmitter", __name__)

class ProgramSubmitter(object):
    '''
    Submits arguments to the cluster or just a general machine. One should
    inherit from this class and implement/override the necessary methods.
    '''

    def __init__(self, argCreator, output_dir, label, timePrefix=False,
            blocking=True, forceExisting=True, numCPU=1, createOutputDir=True,
            ignoreSubmitErrors=False):
        '''
        Create the submitter object.

        Parameters
        ----------
        argCreator : ArgumentCreator
            Argument creator
        output_dir : string
            Output directory to save data to. Unrelated to the actual output
            directory of the program that will run as the submitter cannot know
            how to force the program to output its data to this directory
        appLabel   : string
            Simulation run label. Each simulation will have its own directory
            created that will contain this label
        timePrefix : bool
            Whether to prefix the label with time before the output directory
            will be created.
        forceExisting : bool, optional
            When True, will not raise an error if the output directory already exists.
        ignoreSubmitErrors : bool, optional
            If submission of one process fails, ignore this and continue.
            Otherwise raise an exception.
        '''
        self._ac = argCreator
        self._output_dir = output_dir
        self._label      = label
        if (label == ''):
            self._timePrefix = True
        else:
            self._timePrefix = timePrefix
        self._forceExisting = forceExisting
        self._ignoreSubmitErrors = ignoreSubmitErrors

        self.createOutputDir = createOutputDir
        if self.createOutputDir:
            self._outputDir = self._createOutputDir()
            try:
                os.makedirs(self.outputDir())
            except OSError as e:
                if (not self._forceExisting):
                    raise e
                if (e.errno == errno.EEXIST):
                    log_warn("root.submitters", "Output directory already exists. This"
                        + " might overwrite files.")
                else:
                    raise e
        else:
            self._outputDir = '__not_set__'
            progSLogger.info('Not making output directory. Disabled by the user.')

        # process management
        self._pList = []
        self._blocking = blocking
        self._numCPU = numCPU


    def _createOutputDir(self):
        if (self._timePrefix):
            nowStr = datetime.now().strftime("%Y-%m-%dT%H-%M-%S_")
        else:
            nowStr = ''
        ret = "{0}/{1}{2}".format(self._output_dir, nowStr, self._label)
        return ret.replace(" ", "_")

    def outputDir(self):
        return self._outputDir


    def RunProgram(self, args, job_num):
        '''Should be overridden; args: argument string'''
        raise NotImplementedError()


    def _wait(self):
        '''Remove the first finished process from the process pool.'''
        if (self._blocking == False or len(self._pList) < self._numCPU):
            return

        while True:
            for p_idx, p in enumerate(self._pList):
                errno = p.poll()
                if errno is None:
                    time.sleep(10e-3)  # Sleep 10 ms
                else:
                    self._pList.pop(p_idx)
                    if not self._ignoreSubmitErrors and (errno is not None and errno != 0):
                        raise SubmitError()
                    return


    def _addProcess(self, p):
        if (p is None):
            return
        self._pList.append(p)
        print("  Process list length: {0}".format(len(self._pList)))


    def _cleanup(self):
        print("Waiting for the remaining processes...")
        while (len(self._pList) != 0):
            p = self._pList.pop(0)
            p.wait()

    def waitForProcesses(self):
        self._cleanup()

    def submitAll(self, startJobNum, repeat=1, dry_run=False, filter=None):
        '''
        Submits all the generated jobs. Parameters:
          startJobNum   Start job number index
          repeat        Number of repeats for each parameters dictionary
        '''

        prt = []
        curr_job_num = startJobNum
        for it in range(self._ac.listSize()):
            for rep in range(repeat):
                self._wait()
                if (filter is None) or (filter == it):
                    print "Submitting simulation " + str(it)
                    p = self.RunProgram(self._ac.getArgString(it, curr_job_num ),
                            curr_job_num, dry_run)
                    self._addProcess(p)
                prt.append((curr_job_num, self._ac.getPrintArgString(it)))
                curr_job_num += 1

        # Cleanup
        self._cleanup()

        if (self._ac.printout):
            self._printStr = self.getPrintoutString(prt)
            print self._printStr
        else:
            self._printStr = None



    def saveIterParams(self, iterParams, dimension_labels, dimensions,
                       fileName='iterparams.h5', dry_run=False):
        '''
        Save iterated parameters.

        .. todo::
            Check the consistency of iterParams, dimension_labels and
            dimensions.

        Parameters
        ----------
        iterParams : dict
            A dictionary of iterated parameters. {<name> : np.ndarray}
        dimension_labels : list or tuple of strings
            Labels of the dimensions. These have to reflect ``iterParams``.
        dimensions : list of ints
            Dimension for each label, must be the same length as
            ``dimension_labels``.
        fileName : str, optional
            An output file name, with an extension supported by the data_storage
            package.
        dry_run : bool
            If ``True`` perform only the dry run.
        '''
        if len(dimension_labels) != len(dimensions):
            raise ValueError("len(dimension_labels) != len(dimensions)")

        filePath = os.path.join(self.outputDir(), fileName)
        log_info('root.submitters',
                 'Saving parameter iteration data to: {0}'.format(filePath))
        if dry_run:
            log_info('root.submitters', 'Dry run: not performing the save '
                     'actually.')
        else:
            o = DataStorage.open(filePath, 'w')
            o['iterParams'] = iterParams
            o['dimension_labels'] = list(dimension_labels)
            o['dimensions'] = list(dimensions)
            o.close()

    def _saveAllOptions(self):
        '''Save the iterated options.'''
        iterFileName = self.outputDir() + '/iterparams.h5'
        log_info('root.submitters', 'Saving parameter iteration data: {0}'.format(iterFileName))
        iterOut = DataStorage.open(iterFileName, 'w')
        iterOut['items'] = self._ac.optionList()
        iterOut.close()


    def submitOne(self, it, startJobNum, repeat=1):
        curr_job_num = startJobNum
        for rep in range(repeat):
            print "Submitting simulation " + str(it)
            p = self.RunProgram(self._ac.getArgString(it, curr_job_num),
                    curr_job_num, dry_run)
            self._addProcess(p)
            curr_job_num += 1


    def getPrintoutString(self, prt):
        res = ""
        for vals in prt:
            res += "{0} {1}\n".format(vals[0], vals[1])
        return res


    ## Export parameter sweep information to a text file.
    #
    # The format is:
    #   jobNum param1 value param2 value ...
    #
    # @param fname File name to write to
    #
    def exportIterParams(self, fname):
        f = open(fname, 'w')
        f.write(self._printStr)
        f.close()



class GenericSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a generic multiprocessor machine, either by blocking them, or
    concurrently.
    '''
    def __init__(self, argCreator, appName, output_dir, label, **kw):
        ProgramSubmitter.__init__(self, argCreator, output_dir, label, **kw)
        self._appName  = appName

    def RunProgram(self, args, job_num, dry_run):
        cmdStr = self._appName + ' ' + args
        log_info('root.submitters', 'Submitting command: \n' + cmdStr)
        if dry_run:
            p = None
        else:
            p = subprocess.Popen(cmdStr, shell=True)
        return p



class QsubSubmitter(ProgramSubmitter):
    '''
    Submit jobs on a machine that supports qsub command (cluster)
    '''
    def __init__(self, argCreator, scriptName, qsub_params, output_dir, label,
            **kw):
        ProgramSubmitter.__init__(self, argCreator, output_dir, label, **kw)
        self._scriptName  = scriptName
        self._qsub_params = qsub_params

    def RunProgram(self, args, job_num, dry_run):
        cmdStr = 'qsub ' + self._qsub_params + ' ' + \
                "-o " + self.outputDir() + ' ' + \
                "-N job{0:05} ".format(job_num) + \
                self._scriptName + ' ' + args
        log_info('root.submitters', 'Submitting command: \n' + cmdStr)
        if dry_run:
            p = None
        else:
            p = subprocess.Popen(cmdStr, shell=True)
        return p
