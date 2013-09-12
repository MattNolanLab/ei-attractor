#
#   factory.py
#
#   Submitting factory.
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

from python import WorkstationSubmitter, ClusterSubmitter

class SubmitterFactory(object):
    '''
    Submitting factory.
    '''

    @classmethod
    def getSubmitter(cls, argCreator, appName, envType='guess',
            rtLimit="00:05:00", output_dir=".", label='', timePrefix=True,
            blocking=True, numCPU=1):
        '''
        Get the appropriate factory object, depending on the ``envType```
        parameter. By default, try to guess the best submitter, but this can be
        overriden.

        Parameters
        ----------
        argCreator : ArgumentCreator

        appName : string
            Program to run

        envType : string, optional
            Envirnoment type. Can be one of ``guess``, ``workstation``,
            ``cluster``.

        rtLimit : string
            Run time limit in environments that support it.

        output_dir : string, optional
            Set output directory for the submitting information. Note that this
            does not set the output directory for your script!

        label : string, optional
            Label of the simulation set. Every simulation set will have its own
            directory created for it, every time the submitAll method is
            called. Note that label = '' implies appendTime == True

        appendTime : bool, optional
            Whether to append a time snapshot to the simulatin label. See note
            in label description for more details.

        blocking : bool, optional
            If True, the submitter will wait until the simulation script
            finishes and only then run the next one. Note that this might have
            no effect on systems where the simulation gets detached from the
            submission (e.g. qsub commands)
        '''
        Scls = None
        if (envType == 'guess'):
            raise NotImplementedError()
        elif (envType == 'workstation'):
            Scls =  WorkstationSubmitter
        elif (envType == 'cluster'):
            Scls = ClusterSubmitter
        else:
            raise ValueError("Unsupported envType parameter.")

        return Scls(argCreator, appName, rtLimit, output_dir, label,
                timePrefix, blocking, numCPU)


