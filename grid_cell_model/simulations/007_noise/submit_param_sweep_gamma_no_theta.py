#!/usr/bin/env python
'''
These are the same simulations as submit_param_sweep_gamma, except that the
external theta current is replaced by a current with a constant amplitude that
delivers the same charge per theta cycle.
'''
from param_sweep    import submitParamSweep, SubmissionParser
from default_params import defaultParameters as dp

parser = SubmissionParser()
o = parser.parse_args()

for noise_sigma in parser.noise_sigmas:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = o.env
    simRootDir  = o.where
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = '../common/simulation_stationary.py'
    rtLimit     = '01:30:00' if o.rtLimit is None else o.rtLimit
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

    # No theta parameters
    p['Iext_e_const'] = 487.5   # pA
    p['Iext_i_const'] = 212.5   # pA
    p['Iext_e_theta'] = 0       # pA
    p['Iext_i_theta'] = 0       # pA


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################
    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel, appName,
                     rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
                     rc=parser.rowcol)
