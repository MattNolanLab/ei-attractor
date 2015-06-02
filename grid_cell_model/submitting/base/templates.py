'''Basic simulation templates.'''
from __future__ import absolute_import, print_function, division

from ..noise.templates import SimulationTemplate
from .parsers import GenericSubmissionParser
from ..factory import SubmitterFactory
from ..arguments import ArgumentCreator


class DemoSimulation(SimulationTemplate):
    '''A demonstration of submitting a grid cell network simulation.

    Parameters
    ----------
    app_name : str
        Path to application to run. The actual simulation script/program to
        run.
    sim_label : str
        A lable for the simulation. A directory at the end of the chain of the
        output directories will be created and this will be passed on to the
        simulation script as the ``output_dir`` parameter.
    default_parameters : ConfigObj or dict-like
        A configuration object (could potentially be a dictionary) for default
        parameters. This should be of type ``ConfigObj``.
    user_parameters : dict-like
        Extra parameters that will override the default parameters. These will
        be merged into the ``default_parameters``.
    '''
    def __init__(self, app_name, sim_label, default_parameters, user_parameters=None):
        super(DemoSimulation, self).__init__(app_name,
                                             GenericSubmissionParser(),
                                             default_parameters,
                                             user_parameters)
        self._sim_label = sim_label

    @property
    def sim_label(self):
        '''Simulation label - data should be saved there.'''
        return self._sim_label

    def run(self):
        if self.parser.options is None:
            self.parser.parse_args()
        o = self.parser.options

        p = self._dp.copy()
        p['time']        = 10e3 if o.time is None else o.time  # ms

        p['nthreads']    = o.nthreads
        p['ntrials']     = o.ntrials
        p['verbosity']   = o.verbosity

        p.update(self._up)

        # Submitting
        ac = ArgumentCreator(p, printout=True)

        numRepeat   = 1
        submitter = SubmitterFactory.getSubmitter(ac,
                                                  self._app_name,
                                                  envType=o.env,
                                                  rtLimit=o.rtLimit,
                                                  output_dir=o.where,
                                                  label=self.sim_label,
                                                  blocking=True,
                                                  timePrefix=False,
                                                  numCPU=o.nCPU)
        ac.setOption('output_dir', submitter.outputDir())
        startJobNum = 0
        submitter.submitAll(startJobNum, numRepeat, dry_run=o.dry_run)
