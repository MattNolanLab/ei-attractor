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
from scipy.io           import loadmat, savemat
from matplotlib.pyplot  import *

from tools              import *

import numpy
import numpy as np
import pp
import tools


jobStart = 0
jobEnd   = 10000


def spikeStatistics(spikeTimes, theta_start_t, theta_freq, T):
    '''
    Count
      - number of spikes per each theta cycle and return as an array.
      - phase of first spike with respect to theta
    '''
    theta_T = (1. / theta_freq)
    spikes_offset = spikeTimes - theta_start_t
    spikes_theta_no = numpy.array(spikes_offset // theta_T, dtype=int)
    nThetas = (T - theta_start_t) // theta_T
    spikesPerEachTheta = []
    phaseOfFirstSpike = []
    for theta_it in xrange(nThetas):
        # No. of spikes per theta
        spikesTheta = spikes_theta_no == theta_it
        spikesPerEachTheta.append(numpy.count_nonzero(spikesTheta))

        # Time of the first spike, nans removed
        firstSpike = spikesTheta.nonzero()[0]
        if len(firstSpike) > 0:
            phaseOfFirstSpike.append( (spikes_offset[firstSpike[0]] -
                    theta_it*theta_T) / theta_T * 2 * numpy.pi )

    return numpy.mean(spikesPerEachTheta), numpy.mean(phaseOfFirstSpike)
        


def currentStatistics(current, dt, theta_start_t, theta_freq):
    '''
    Remove the portion of the signal before theta_start_t, and then compute:
      * amplitude of the current (simply max) per each theta cycle
      * total charge per theta cycle
    '''
    curr_phase = tools.splitSigToThetaCycles(current[theta_start_t/dt:],
            1./theta_freq, dt)

    curr_max = numpy.max(curr_phase, axis=1)
    curr_min = numpy.min(curr_phase, axis=1)
    charge = tools.getChargeTheta(curr_phase, dt)

    return numpy.mean(curr_max), numpy.mean(curr_min), numpy.mean(charge)


def upSample(sig, dt):
    return numpy.repeat(sig, 2), dt/2.0 


def loadAndDoStats(job_num):
    print 'job_num', job_num

    output_dir = 'output_local/'
    fileNamePrefix = ''
    trial_it = 0
    
    input_fname = "{0}/{1}job{2:05}_trial{3:04}_output.mat".format(output_dir,
            fileNamePrefix, job_num, trial_it)
    inData = scipy.io.loadmat(input_fname, squeeze_me=False)
    
    
    stateMon_times = inData['stateMon_times'][:, 0]
    mon_dt = stateMon_times[1] - stateMon_times[0]
    theta_start_t = float(inData['theta_start_mon_t']) / 1e3
    theta_freq = int(inData['options']['theta_freq'])
    sim_T = float(inData['options']['time']) / 1e3
    
    ################################################################################
    #                              E/I cell statistics
    ################################################################################
    e_spikesPerTheta    = []
    e_phaseOfFirstSpike = []
    e_curr_max = []
    e_curr_min = []
    e_charge   = []

    i_spikesPerTheta    = []
    i_phaseOfFirstSpike = []
    i_curr_max = []
    i_curr_min = []
    i_charge   = []

    for EI_id in [0, 1]:
        if (EI_id == 0):
            # E cells
            spikes          = inData['spikeCell_e']
            Iclamp          = inData['stateMon_Iclamp_e_values']
            Ncells          = len(spikes)
            Ncells_state    = len(Iclamp)

            spikesPerTheta      = e_spikesPerTheta
            phaseOfFirstSpike   = e_phaseOfFirstSpike
            curr_max            = e_curr_max
            curr_min            = e_curr_min
            charge              = e_charge
        elif (EI_id == 1):
            # I cells
            spikes          = inData['spikeCell_i']
            Iclamp          = inData['stateMon_Iclamp_i_values']
            Ncells          = len(spikes)
            Ncells_state    = len(Iclamp)

            spikesPerTheta      = i_spikesPerTheta
            phaseOfFirstSpike   = i_phaseOfFirstSpike
            curr_max            = i_curr_max
            curr_min            = i_curr_min
            charge              = i_charge
        else:
            raise IndexError("EI_id is not 0(E) or 1(I)!")

        for n_it in xrange(Ncells):
            sp, ph = spikeStatistics(spikes[n_it, 0].flatten(), theta_start_t,
                    theta_freq, sim_T)
            spikesPerTheta.append(sp)
            phaseOfFirstSpike.append(ph)


        for n_it in xrange(Ncells_state):
            data, mon_dt_up = upSample(Iclamp[n_it, :],
                    mon_dt)
            c_max, c_min, Q = currentStatistics(data, mon_dt_up, theta_start_t, theta_freq)
            curr_max.append(c_max)
            curr_min.append(c_min)
            charge.append(Q)

    #import pdb; pdb.set_trace()

    return e_spikesPerTheta, \
           e_phaseOfFirstSpike, \
           e_curr_max, \
           e_curr_min, \
           e_charge, \
           i_spikesPerTheta, \
           i_phaseOfFirstSpike, \
           i_curr_max, \
           i_curr_min, \
           i_charge





#Ndim = 10
#spikesPerTheta_all      = numpy.ndarray((Ndim, Ndim, Ndim, Ndim), dtype=object)
#phaseOfFirstSpike_all   = numpy.ndarray((Ndim, Ndim, Ndim, Ndim), dtype=object)

parallel = True

all_data = []

if parallel:
    job_server = pp.Server(ppservers=('jupiter1', 'localhost', 'jupiter3'))
    
    jobs = []
    for job_num in xrange(jobStart, jobEnd):
        jobs.append(job_server.submit(
            loadAndDoStats,
            (job_num,),
            (spikeStatistics, currentStatistics, upSample),
            ("numpy", "scipy.io", "tools")))
    
    for job_num in xrange(jobStart, jobEnd):
        all_data.append(jobs[job_num - jobStart]())
else:
    for job_num in xrange(jobStart, jobEnd):
        all_data.append(loadAndDoStats(job_num))

#job_num = jobStart
#for theta_depth in xrange(Ndim):
#    for Iext_e_amp in xrange(Ndim):
#        for E_coupling in xrange(Ndim):
#            for I_coupling in xrange(Ndim):
#                #if job_num > 0: continue
#
#                spikesPerTheta_all[theta_depth, Iext_e_amp, E_coupling,
#                        I_coupling] = numpy.array(spikesPerTheta)
#                phaseOfFirstSpike_all[theta_depth, Iext_e_amp, E_coupling,
#                        I_coupling] = numpy.array(phaseOfFirstSpike)
#                
#                
#                job_num += 1

all_data = np.array(all_data)

outData = {}
outData['all_data'] = all_data
#outData['spikesPerTheta']       = all_data[:, 0]
#outData['phaseOfFirstSpike']    = all_data[:, 1]
#outData['curr_max']             = all_data[:, 2]
#outData['curr_min']             = all_data[:, 3]
#outData['charge']               = all_data[:, 4]


outDataDir = 'output_local'
savemat(outDataDir + '/param_sweep_data_analysis.mat', outData)

