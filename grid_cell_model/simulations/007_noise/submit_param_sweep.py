#!/usr/bin/env python
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
from default_params         import defaultParameters as p
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator


# Submitting
ENV         = 'workstation'
simRootDir  =  'output_local'
simLabel    =  'tmp_param_sweep'
appName     = 'simulation_stationary.py'
rtLimit     = '00:05:00'
blocking    = True
timePrefix  = False
numRepeat   = 1

#p['time']              = 1199.9e3  # ms
p['time']              = 2e3  # ms
p['nthreads']          = 8
p['ntrials']           = 2

p['noise_sigma']       = 150.0     # pA


###############################################################################
ac = ArgumentCreator(p, printout=True)

## Range of parameters around default values
## Let's choose a 10% jitter around the default values
#Ndim        = 10     # Number of values for each dimension
#jitter_frac = 0.1    # Fraction
#Iext_e_amp_default = p['Iext_e_const'] + p['Iext_e_theta']
#
#jitter_frac_arr = np.linspace(1.0-jitter_frac, 1.0+jitter_frac, Ndim)
#
#theta_depth_range = jitter_frac_arr*p['Iext_e_const']
#Iext_e_amp_range  = jitter_frac_arr*Iext_e_amp_default
#
#Iext_e_const_arr     = []
#Iext_e_theta_arr     = []
#g_AMPA_total_arr     = []
#g_GABA_total_arr     = []
#g_uni_GABA_total_arr = []
#for theta_depth in theta_depth_range:
#    for Iext_e_amp in Iext_e_amp_range:
#        for E_coupling in jitter_frac_arr:
#            for I_coupling in jitter_frac_arr:
#                Iext_e_const_arr.append(theta_depth)
#                Iext_e_theta_arr.append(Iext_e_amp - theta_depth)
#                g_AMPA_total_arr.append(E_coupling*p['g_AMPA_total'])
#                g_GABA_total_arr.append(I_coupling*p['g_GABA_total'])
#                g_uni_GABA_total_arr.append(I_coupling*p['g_uni_GABA_total'])
#
#
#iterparams = {
#        'Iext_e_const'      : Iext_e_const_arr,
#        'Iext_e_theta'      : Iext_e_theta_arr,
#        'g_AMPA_total'      : g_AMPA_total_arr,
#        'g_GABA_total'      : g_GABA_total_arr,
#        'g_uni_GABA_total'  : g_uni_GABA_total_arr
#}
#ac.insertDict(iterparams, mult=False)


###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix);
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
