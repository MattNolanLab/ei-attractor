#
#   submit_param_sweep.py
#
#   Submit job(s) to the cluster/workstation: parameter sweep (noise)
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

from default_params import defaultParameters
from common         import *
from scipy.io       import savemat

import logging  as lg
import numpy    as np


lg.basicConfig(level=lg.DEBUG)


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['output_dir'] = 'output_local'

parameters['time']              = 6e3   # ms
parameters['theta_start_mon_t'] = 1.0e3   # ms

parameters['prefDirC_e']        = 0
parameters['prefDirC_i']        = 0

parameters['theta_noise_sigma'] = 0         # pA


parameters['noise_sigma']       = 0.0       # mV



startJobNum = 1
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_param_sweep.py'
blocking            = True

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_param_sweep.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:10:00"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

# Range of parameters around default values
# Let's choose a 10% jitter around the default values
Ndim        = 2     # Number of values for each dimension
jitter_frac = 0.1    # Fraction
Iext_e_amp_default = parameters['Iext_e_const'] + parameters['Iext_e_theta']

jitter_frac_arr = np.linspace(1.0-jitter_frac, 1.0+jitter_frac, Ndim)

theta_depth_range = jitter_frac_arr*parameters['Iext_e_const']
Iext_e_amp_range  = jitter_frac_arr*Iext_e_amp_default

Iext_e_const_arr     = []
Iext_e_theta_arr     = []
g_AMPA_total_arr     = []
g_GABA_total_arr     = []
g_uni_GABA_total_arr = []
for theta_depth in theta_depth_range:
    for Iext_e_amp in Iext_e_amp_range:
        for E_coupling in jitter_frac_arr:
            for I_coupling in jitter_frac_arr:
                Iext_e_const_arr.append(theta_depth)
                Iext_e_theta_arr.append(Iext_e_amp - theta_depth)
                g_AMPA_total_arr.append(E_coupling*parameters['g_AMPA_total'])
                g_GABA_total_arr.append(I_coupling*parameters['g_GABA_total'])
                g_uni_GABA_total_arr.append(I_coupling*parameters['g_uni_GABA_total'])


iterparams = {
        'Iext_e_const'      : Iext_e_const_arr,
        'Iext_e_theta'      : Iext_e_theta_arr,
        'g_AMPA_total'      : g_AMPA_total_arr,
        'g_GABA_total'      : g_GABA_total_arr,
        'g_uni_GABA_total'  : g_uni_GABA_total_arr
}
ac.insertDict(iterparams, mult=False, printout=True)


if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=True)

# Export the iterparams
savemat(parameters['output_dir'] + '/param_sweep_iterparams.mat', iterparams, oned_as='row')

