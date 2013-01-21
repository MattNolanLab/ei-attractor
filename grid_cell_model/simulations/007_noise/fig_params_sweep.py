#
#   fig_param_sweep.py
#
#   Parameter sweep: data analysis and figures
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
from scipy              import linspace
from scipy.io           import loadmat, savemat
from matplotlib.pyplot  import *

import numpy as np


def countSpikesPerTheta(spikeTimes, theta_start_t, theta_freq):
    '''
    Count number of spikes per each theta cycle and return as an array.
    '''
    spikes_offset = spikeTimes - theta_start_t
    spikes_theta_no = np.array(spikes_offset // (1. / theta_freq), dtype=int)
    nThetas = int(np.max(spikes_theta_no))
    spikesPerEachTheta = []
    for theta_it in xrange(nThetas):
        spikesPerEachTheta.append(np.count_nonzero(spikes_theta_no == theta_it))

    return np.array(spikesPerEachTheta)
        
    


output_dir = 'output'
fileNamePrefix = ''
trial_it = 0
jobStart = 0

Ndim = 10
spikesPerTheta_all = np.ndarray((Ndim, Ndim, Ndim, Ndim))


job_num = jobStart
for theta_depth in xrange(Ndim):
    for Iext_e_amp in xrange(Ndim):
        for E_coupling in xrange(Ndim):
            for I_coupling in xrange(Ndim):
                print 'job_num', job_num

                input_fname = "{0}/{1}job{2:05}_trial{3:04}_output.mat".format(output_dir,
                        fileNamePrefix, job_num, trial_it)
                inData = loadmat(input_fname, squeeze_me=True)
                
                
                stateMon_times = inData['stateMon_times']
                mon_dt = stateMon_times[1] - stateMon_times[0]
                theta_start_t = int(inData['theta_start_mon_t']) / 1e3
                theta_freq = int(inData['options']['theta_freq'])
                
                ################################################################################
                #                      Number of spikes per theta cycle
                ################################################################################
                Ncells_e = len(inData['spikeCell_e'])
                spikesPerTheta = []
                for n_it in xrange(Ncells_e):
                    spikesPerTheta.append(countSpikesPerTheta(inData['spikeCell_e'][n_it],
                        theta_start_t, theta_freq))
                spikesPerTheta = np.array(spikesPerTheta)

                spikesPerTheta_all[theta_depth, Iext_e_amp, E_coupling,
                        I_coupling] = np.mean(np.mean(spikesPerTheta, axis=1))
                
                
                job_num += 1


outData = {}
outData['spikesPerTheta_all'] = spikesPerTheta_all
savemat(output_dir + '/param_sweep_data_analysis.mat', outData)

