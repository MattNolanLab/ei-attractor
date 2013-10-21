#
#   fig_gamma_vel_freq.py
#
#   Frequency of gamma dependent on running speed?
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
import numpy as np

from scipy.io           import loadmat
from matplotlib.pyplot  import *
from tables             import *

from analysis.grid_cells import extractSpikePositions2D
from analysis.signal     import butterBandPass, CWT


jobRange = [100, 100]
trialNum = 0
dumpNum = 4

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


f_start = 40
f_stop   = 200

pA = 1e-12


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"


for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    print 'jobNum: ' + str(jobNum)

    fileName = fileNameTemp +  '_output.mat'
    fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                                trialNum, dumpNum)
    try:
        data = loadmat(fileName)
    except:
        print "warning: could not open: " + fileName
        continue

    pos_x  = data['pos_x'].ravel()
    pos_y  = data['pos_y'].ravel()
    rat_dt = data['dt'][0][0]
    velocityStart = data['velocityStart'][0][0]

    times    = data['stateMon_Iclamp_e_times'].ravel()
    sim_dt   = times[1] - times[0]
    Iclamp_e = data['stateMon_Iclamp_e_values'][:, 0].ravel()/pA

    vel_x = np.diff(pos_x)/rat_dt
    vel_y = np.diff(pos_y)/rat_dt
    spd = np.sqrt(vel_x**2 + vel_y**2)

    Iclamp_e_bpass = butterBandPass(Iclamp_e, sim_dt, f_start, f_stop)
    figure()
    plot(times, Iclamp_e_bpass/pA)

    #NFFT = 256
    #noverlap = NFFT/4
    #Fs   = 1./sim_dt
    #figure()
    #specgram(Iclamp_e_bpass, NFFT, Fs, noverlap=noverlap)
    #ylim([0, f_stop])
    #xlim([0, 5])

    figure()
    subplot(2, 1, 1)
    ts = 280
    te = 299
    t_start = ts - times[0]
    t_end   = te - times[0]
    time_slice = times[t_start/sim_dt:t_end/sim_dt]
    Iclamp_w, Iclamp_f = CWT(Iclamp_e_bpass[t_start/sim_dt:t_end/sim_dt], sim_dt, f_stop)
    T, F = np.meshgrid(time_slice, Iclamp_f)
    pcolormesh(T, F, Iclamp_w, edgecolors='None', cmap=get_cmap('jet'))
    minorticks_on()
    ylabel('F (Hz)')
    xlim([ts, te])


    # These must be shifted to accomodate for the shift in velocity input start time
    t_vel_start = ts - velocityStart*1e-3
    t_vel_end   = te - velocityStart*1e-3
    spd_slice = spd[t_vel_start/rat_dt: t_vel_end/rat_dt]
    subplot(2, 1, 2)
    plot(np.arange(ts, te, rat_dt), spd_slice)
    xlim([ts, te])
    ylim([0, 100])
    xlabel('time (s)')
    ylabel('rat speed (cm/s)')


    # Log scale
    Iclamp_w = 10 * np.log10(Iclamp_w)
    figure()
    subplot(2, 1, 1)
    pcolormesh(T, F, Iclamp_w, edgecolors='None', cmap=get_cmap('jet'))
    minorticks_on()
    ylabel('F (Hz)')
    xlim([ts, te])


    # These must be shifted to accomodate for the shift in velocity input start time
    subplot(2, 1, 2)
    plot(np.arange(ts, te, rat_dt), spd_slice)
    xlim([ts, te])
    ylim([0, 100])
    xlabel('time (s)')
    ylabel('rat speed (cm/s)')



    show()

    #h5file = openFile(fileName + '_' + spikeType + '_igor.h5', mode = "w", title =
    #        "Bump velocity export figures")
    #
    #neuronPos_x, neuronPos_y, max_i = extractSpikePositions2D(spikes, pos_x, pos_y, rat_dt)
    #h5file.createArray(h5file.root, 'rat_pos_x', pos_x[0:max_i+1])
    #h5file.createArray(h5file.root, 'rat_pos_y', pos_y[0:max_i+1])
    #h5file.createArray(h5file.root, 'neuronPos_x', neuronPos_x)
    #h5file.createArray(h5file.root, 'neuronPos_y', neuronPos_y)

    #rateMap[rateMap.mask == 1] = np.nan
    #h5file.createArray(h5file.root, 'rateMap', rateMap)
    #corr[corr.mask==1] = np.nan
    #h5file.createArray(h5file.root, 'corrMap', corr)

    #h5file.close()

