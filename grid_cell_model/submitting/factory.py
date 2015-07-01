'''Submitting factory.

Determines which simulation submitter to use, depending on the environment
settings.
'''
from .python import WorkstationSubmitter, ClusterSubmitter

class SubmitterFactory(object):
    '''
    Submitting factory.
    '''
    @staticmethod
    def getSubmitter(argCreator, appName, envType='guess', output_dir=".",
            label='', **kw):
        '''
        Get the appropriate factory object, depending on the ``envType```
        parameter. By default, try to guess the best submitter, but this can be
        overriden.

        Parameters
        ----------
        argCreator : ArgumentCreator
            Argument creator to use with the submitter.
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
            called. Note that label = '' implies appendTime == True>
        appendTime : bool, optional
            Whether to append a time snapshot to the simulatin label. See note
            in label description for more details.
        blocking : bool, optional
            If True, the submitter will wait until the simulation script
            finishes and only then run the next one. Note that this might have
            no effect on systems where the simulation gets detached from the
            submission (e.g. qsub commands).
        '''
        Scls = None
        if (envType == 'guess'):
            raise NotImplementedError()
        elif (envType == 'workstation'):
            Scls =  WorkstationSubmitter
            kw.pop('rtLimit', None)            # ignored here
            kw.pop('extra_qsub_params', None)  # ignored here
        elif (envType == 'cluster'):
            Scls = ClusterSubmitter
        else:
            raise ValueError("Unsupported envType parameter.")

        return Scls(argCreator, appName, output_dir, label, **kw)


