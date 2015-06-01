#!/usr/bin/env python
'''Parameter exploration simulations in which one parameter is being changed.

Run a network for a short time and save data. These simulations are meant to be
run for interactive tuning of the network.
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.noise import SingleParameterSweepParser
from default_params import defaultParameters as dp

parser = SingleParameterSweepParser()
parser.add_argument('--Ivel', type=float,
                    help='Velocity input (pA). Default is 50 pA.')
o = parser.parse_args()

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
    p['nthreads']    = 1
    p['ntrials']     = o.ntrials
    p['verbosity']   = o.verbosity
    p['Ivel']        = 50. if o.Ivel is None else o.Ivel  # mA

    p['g_AMPA_total'] = 290.    # nS
    p['g_GABA_total'] = 410.    # nS
    p['g_EE_total']   = 500.    # nS

    # No theta parameters
    p['Iext_e_const'] = 400.0   # pA
    p['Iext_i_const'] = 175.0   # pA
    p['Iext_e_theta'] = 0       # pA
    p['Iext_i_theta'] = 0       # pA

    # Here, no PC inputs
    p['pc_start_max_rate'] = .0  # Hz

    p['prefDirC_e'] = 0.

    explored_param = np.arange(o.param_start, o.param_stop + o.param_step,
                               o.param_step)
    dummy = explored_param * 0. + p['prefDirC_e']
    iterparams = {
        'prefDirC_e': dummy,
        o.explored_param: explored_param
    }
    dimension_labels = ['prefDirC_e', o.explored_param]
    dimensions = [1, len(explored_param)]
    ac = ArgumentCreator(p, printout=True)
    ac.insertDict(iterparams, mult=False)

    ###########################################################################
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix,
        numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
    submitter.saveIterParams(iterparams, dimension_labels, dimensions,
                             dry_run=dry_run)
