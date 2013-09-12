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

from analysis.signal    import splitSigToThetaCycles, getChargeTheta
from analysis.spikes    import firingRate, multipleFiringRate

import numpy
import numpy as np
import pp


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
        spikesTheta = (spikes_theta_no == theta_it)
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
    curr_phase = splitSigToThetaCycles(current[theta_start_t/dt:],
            1./theta_freq, dt)

    curr_max = numpy.max(curr_phase, axis=1)
    curr_min = numpy.min(curr_phase, axis=1)
    charge = getChargeTheta(curr_phase, dt)

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


def deleteNans(arr):
    '''
    Removes all nan values from a numpy array
    '''
    return np.delete(arr, np.nonzero(np.isnan(arr)))


def getFiringRate(jobNum):
    '''
    Load one job with number jobNum and process the E and I population firing
    rate. Return as tuple (e, i)
    '''
    output_dir = 'output_local/'
    fileNamePrefix = ''
    trial_it = 0
    input_fname = "{0}/{1}job{2:05}_trial{3:04}_output.mat".format(output_dir,
            fileNamePrefix, jobNum, trial_it)

    inData = scipy.io.loadmat(input_fname, squeeze_me=False, variable_names=
            ["spikeCell_e", "spikeCell_i", "options"])

    thetaStartTime = int(inData['options']['theta_start_mon_t']) / 1e3
    time = int(inData['options']['time']) / 1e3
    winLen = 1  # sec
    y_dim = np.sqrt(3) / 2.0
    Ne_x = int(inData['options']['Ne'])
    Ne_y = int(np.ceil(Ne_x * y_dim)) // 2 * 2
    Ni_x = int(inData['options']['Ni'])
    Ni_y = int(np.ceil(Ni_x * y_dim)) // 2 * 2

    print Ne_x, Ne_y, Ni_x, Ni_y

    firingRate_e = multipleFiringRate(inData['spikeCell_e'].flatten(), time - winLen, time)
    firingRate_i = multipleFiringRate(inData['spikeCell_i'].flatten(), time - winLen, time)

    return (np.reshape(firingRate_e, (Ne_y, Ne_x)), np.reshape(firingRate_i,
        (Ni_y, Ni_x)))
    


################################################################################
#                                Main run
################################################################################

Ndim = 10


preprocess  = False
parallel    = True

outDataDir = 'output_local'


if preprocess == True:
    # Load all the simulation data one by one and process spikes and currents
    print "Preprocessing data to produce aggregate information"

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

    all_data = np.array(all_data)
    
    outData = {}
    outData['all_data'] = all_data
    
    
    savemat(outDataDir + '/param_sweep_data_analysis.mat', outData)


else:
    # Load data analysis data and do plots/exports
    input_fname = outDataDir + '/param_sweep_data_analysis.mat'
    iterparams_fname = outDataDir + '/param_sweep_iterparams.mat'
    print "Loading \'" + input_fname + "\' and plotting/exporting plots"

    inputDict = loadmat(input_fname)
    iterparams = loadmat(iterparams_fname)

    multiDimData     = numpy.ndarray((Ndim, Ndim, Ndim, Ndim), dtype=object)
    
    job_num = jobStart
    for theta_depth in xrange(Ndim):
        for Iext_e_amp in xrange(Ndim):
            for E_coupling in xrange(Ndim):
                for I_coupling in xrange(Ndim):
                    multiDimData[theta_depth, Iext_e_amp, E_coupling,
                            I_coupling] = tuple(inputDict['all_data'][job_num, :])
                    job_num += 1

    theta_depth_arr = np.reshape(iterparams['Iext_e_const'].flatten(),
            (Ndim, Ndim, Ndim, Ndim))
    


    from tables import *
    h5file = openFile('output/param_sweep_igor_export.h5', mode='w', title = 
            "Parameter sweep export")


    ###########################################################################
    #                             Theta-depth
    ###########################################################################
    theta_depth_gr = h5file.createGroup("/", 'theta_depth', 'Change of theta depth')

    theta_depth_x = theta_depth_arr[:, 0.5*Ndim, 0.5*Ndim, 0.5*Ndim]

    example_it = int(Ndim / 2)

    e_spikesPerTheta         = []
    e_spikesPerTheta_std     = []
    e_phaseOfFirstSpike      = []
    e_phaseOfFirstSpike_std  = []
    e_curr_max               = []
    e_curr_max_std           = []
    i_spikesPerTheta         = []
    i_spikesPerTheta_std     = []
    i_curr_min               = []
    i_curr_min_std           = []

    for theta_depth_it in xrange(Ndim):
        e_spikes = multiDimData[theta_depth_it, 0.5*Ndim, 0.5*Ndim,
                0.5*Ndim][0].flatten()
        e_spikesPerTheta.append(np.mean(e_spikes))
        e_spikesPerTheta_std.append(np.std(e_spikes))

        e_phase_nonan = deleteNans(multiDimData[theta_depth_it, 0.5*Ndim, 0.5*Ndim,
                0.5*Ndim][1].flatten())
        e_phaseOfFirstSpike.append(np.mean(e_phase_nonan))
        e_phaseOfFirstSpike_std.append(np.std(e_phase_nonan))

        e_curr = multiDimData[theta_depth_it, 0.5*Ndim, 0.5*Ndim,
                0.5*Ndim][2].flatten()
        e_curr_max.append(np.mean(e_curr))
        e_curr_max_std.append(np.std(e_curr))

        i_spikes = multiDimData[theta_depth_it, 0.5*Ndim, 0.5*Ndim,
                0.5*Ndim][5].flatten()
        i_spikesPerTheta.append(np.mean(i_spikes))
        i_spikesPerTheta_std.append(np.std(i_spikes))

        i_curr = multiDimData[theta_depth_it, 0.5*Ndim, 0.5*Ndim,
                0.5*Ndim][8].flatten()
        i_curr_min.append(np.mean(i_curr))
        i_curr_min_std.append(np.std(i_curr))

        # Choose examples of data
        if theta_depth_it == example_it:
            e_phaseOfFirstSpike_ex  = e_phase_nonan
            e_spikesPerTheta_ex     = e_spikes
            i_spikesPerTheta_ex     = i_spikes
            e_curr_max_ex           = e_curr
            i_curr_min_ex           = i_curr
            theta_depth_ex          = theta_depth_x[theta_depth_it]

            exJobNum = example_it*10**3 + int(0.5*Ndim)*10**2 + int(0.5*Ndim)*10 +int(0.5*Ndim)
            print "Theta depth example job number is " + str(exJobNum)
            e_firingRate_ex, i_firingRate_ex = getFiringRate(exJobNum)
            startJobNum = 0*10**3 + int(0.5*Ndim)*10**2 + int(0.5*Ndim)*10 +int(0.5*Ndim)
            endJobNum = (Ndim-1)*10**3 + int(0.5*Ndim)*10**2 + int(0.5*Ndim)*10 +int(0.5*Ndim)
            e_firingRate_start, i_firingRate_start = getFiringRate(startJobNum)
            e_firingRate_end, i_firingRate_end = getFiringRate(endJobNum)

    
    e_spikesPerTheta_sem = np.array(e_spikesPerTheta_std) / np.sqrt(Ndim)
    e_phaseOfFirstSpike_sem = np.array(e_phaseOfFirstSpike_std) / np.sqrt(Ndim)
    i_spikesPerTheta_sem = np.array(i_spikesPerTheta_std) / np.sqrt(Ndim)


    errorbar(theta_depth_x, e_spikesPerTheta, yerr=e_spikesPerTheta_sem, fmt='o')
    hold('on')
    errorbar(theta_depth_x, i_spikesPerTheta, yerr=i_spikesPerTheta_sem, fmt='o')
#    ylim([0, 2])


    h5file.createArray(theta_depth_gr, "e_spikesPerTheta",     e_spikesPerTheta)
    h5file.createArray(theta_depth_gr, "e_spikesPerTheta_std", e_spikesPerTheta_std)
    h5file.createArray(theta_depth_gr, "e_spikesPerTheta_sem", e_spikesPerTheta_sem)

    h5file.createArray(theta_depth_gr, "e_phaseOfFirstSpike",     e_phaseOfFirstSpike)
    h5file.createArray(theta_depth_gr, "e_phaseOfFirstSpike_std", e_phaseOfFirstSpike_std)
    h5file.createArray(theta_depth_gr, "e_phaseOfFirstSpike_sem", e_phaseOfFirstSpike_sem)

    h5file.createArray(theta_depth_gr, "e_curr_max",     e_curr_max)
    h5file.createArray(theta_depth_gr, "e_curr_max_std", e_curr_max_std)

    h5file.createArray(theta_depth_gr, "i_spikesPerTheta",     i_spikesPerTheta)
    h5file.createArray(theta_depth_gr, "i_spikesPerTheta_std", i_spikesPerTheta_std)
    h5file.createArray(theta_depth_gr, "i_spikesPerTheta_sem", i_spikesPerTheta_sem)

    h5file.createArray(theta_depth_gr, "i_curr_min",     i_curr_min)
    h5file.createArray(theta_depth_gr, "i_curr_min_std", i_curr_min_std)

    h5file.createArray(theta_depth_gr, "theta_depth", theta_depth_x)

    h5file.createArray(theta_depth_gr, "e_spikesPerTheta_ex",    e_spikesPerTheta_ex)
    h5file.createArray(theta_depth_gr, "e_phaseOfFirstSpike_ex", e_phaseOfFirstSpike_ex)
    h5file.createArray(theta_depth_gr, "e_firingRate_ex",        e_firingRate_ex)
    h5file.createArray(theta_depth_gr, "e_firingRate_start",     e_firingRate_start)
    h5file.createArray(theta_depth_gr, "e_firingRate_end",       e_firingRate_end)

    h5file.createArray(theta_depth_gr, "e_curr_max_ex",          e_curr_max_ex)
    h5file.createArray(theta_depth_gr, "i_spikesPerTheta_ex",    i_spikesPerTheta_ex)
    h5file.createArray(theta_depth_gr, "i_firingRate_ex",        i_firingRate_ex)
    h5file.createArray(theta_depth_gr, "i_curr_min_ex",          i_curr_min_ex)
    h5file.createArray(theta_depth_gr, "theta_depth_ex",         theta_depth_ex)
    h5file.createArray(theta_depth_gr, "i_firingRate_start",     i_firingRate_start)
    h5file.createArray(theta_depth_gr, "i_firingRate_end",       i_firingRate_end)




#show()


