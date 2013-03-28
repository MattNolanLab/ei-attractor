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
            rtLimit="00:05:00", output_dir=".", blocking=True):
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

        return Scls(argCreator, appName, rtLimit, output_dir, blocking)


