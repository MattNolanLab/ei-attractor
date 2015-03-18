#!/usr/bin/env python
'''Run a network for a short time; save data and plot activity.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator

from param_sweep import SubmissionParserBase
from default_params import defaultParameters as dp

parser = SubmissionParserBase()
parser.add_argument('--nthreads', type=int, default=1,
                    help='Number of simulation threads.')
parser.add_argument('--Ivel', type=float,
                    help='Velocity input (pA). Default is 50 pA.')
parser.add_argument('--ee_connections', type=int, choices=[0, 1], default=0,
                    help=('Whether to use E-->E connections (this makes E-->I '
                          'and I-->E flat'))
o = parser.parse_args()

def set_up_connection_params(opts, params):
    '''Set up parameters for the synaptic profiles.

    Parameters
    ----------
    o : Struct-like
        Command line options.
    p : dict-like
        Network parameters.

    Returns
    -------
    p : dict-like
        In-situ modified ``p``.
    '''
    if opts.ee_connections:
        params['EI_flat'] = 1
        params['IE_flat'] = 1
        params['use_EE']  = 1
        params['g_EE_total'] = 4000. #nS
    return p

for noise_sigma in parser.noise_sigmas:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma  # pA

    # Submitting
    ENV         = o.env
    simRootDir  = o.where
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = '../common/simulation_test_network.py'
    rtLimit     = o.rtLimit
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run

    p['master_seed'] = 123456
    p['time']        = 10e3 if o.time is None else o.time  # ms
    p['nthreads']    = o.nthreads
    p['ntrials']     = o.ntrials
    p['verbosity']   = o.verbosity
    p['Ivel']        = 50. if o.Ivel is None else o.Ivel  # mA

    p['g_AMPA_total'] = 180.    # nS
    p['g_GABA_total'] = 200.   # nS

    # No theta parameters
    p['Iext_e_const'] = 400.0   # pA
    p['Iext_i_const'] = 175.0   # pA
    p['Iext_e_theta'] = 0       # pA
    p['Iext_i_theta'] = 0       # pA

    # Here, no PC inputs
    #p['pc_start_max_rate'] = .0  # Hz

    p['prefDirC_e'] = 0.


    p = set_up_connection_params(o, p)

    ###########################################################################
    ac = ArgumentCreator(p, printout=True)
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix,
        numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
