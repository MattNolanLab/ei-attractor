#
#   python.py
#
#   Python submitters
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

from submitters import GenericSubmitter, QsubSubmitter


class ClusterSubmitter(QsubSubmitter):
    '''
    Submit python jobs on a cluster.

    A file named ``cluster_submit.sh`` must be present in the current working
    directory. This script will be run by the qsub command.
    '''

    def __init__(self, argCreator, appName, rtLimit, outputDir,
            label, timePrefix, blocking, numCPU):
        '''
        Initialize.

        Parameters
        ----------
        argCreator : ArgumentCreator object

        appName : string
            A python script that the cluster_submit.sh script should run.

        outputDir : string, optional, default: "."
            Output directory, should be the same output directory of data
            generated by the script.

        rtLimit : string, default: "00:05:00"
            Hard runtime limit, (passed as -l h_rt to qsub)

        blocking : bool
            Whether a running instance should block. Note that this will be
            ignored on a cluster. The commands are always **non-blocking**.

        numCPU : int
            Number of cores to use for OpenMP programs
        '''
        self.submitScript = 'cluster_submit.sh '
        self.default_qsub_params = '-notify ' 
        self.appName = appName
        scriptName = self.submitScript + appName
        self.rtLimit  = rtLimit
        qsub_params = self.default_qsub_params + '-l h_rt=\'' + rtLimit + '\''
        if numCPU > 1:
            qsub_params += ' -pe OpenMP {0}'.format(numCPU)
        elif (numCPU < 1):
            raise ValueError("numCPU must be >= 1.")
        QsubSubmitter.__init__(self, argCreator, scriptName, qsub_params,
                outputDir, label, timePrefix, blocking)




class WorkstationSubmitter(GenericSubmitter):
    '''
    Submit jobs on a workstation.

    This will submit the python command on a generic processor machine.
    '''

    def __init__(self, argCreator, appName, rtLimit, outputDir, label, **kw):
        '''
        Initialize the submitter.

        Parameters
        ----------
        argCreator : ArgumentCreator

        appName : string
            Path to the application that should be run by Python.

        outputDir : string, optional, default: "."
            Output directory, should be the same output directory of data
            generated by the script. Note that currently this does not have any
            effect.

        rtLimit : string
            Hard runtime limit, ignored here.

        blocking : bool
            Whether a running instance should block.

        numCPU : int
            Number of cores to use for OpenMP programs. This parameter is here
            only for interface compatibility. It is ignored in this class
        '''
        self.progName = 'python '
        commandStr = self.progName + appName
        GenericSubmitter.__init__(self, argCreator, commandStr, outputDir,
                label, **kw)

