#
#   tools.py
#
#   Some useful data analysis functions.
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
from scipy.signal import *
from Wavelets import Morlet

from os import system


def butterHighPass(sig, dt, f_pass):
    nyq_f = 1./dt/2
    norm_f_pass = f_pass/nyq_f

    # Low pass filter
    b, a = butter(3, norm_f_pass, btype='high')
    return filtfilt(b, a, sig)


def butterBandPass(sig, dt, f_start, f_stop):
    '''Band pass filter a signal, with f_start and f_stop frequencies'''
    nyq_f = 1./dt/2
    norm_f_start = f_start/ nyq_f
    norm_f_stop  = f_stop / nyq_f
    b, a = butter(3, [norm_f_start, norm_f_stop], btype='band')
    return filtfilt(b, a, sig)


def spikePhaseTrialRaster(spikeTimes, f, start_t=0):
    '''Here assuming that phase(t=0) = 0'''
    spikeTimes -= start_t
    trials = np.floor(f*spikeTimes)
    phases = np.mod(2*np.pi*f*spikeTimes, 2*np.pi)
    times  = np.mod(spikeTimes, 1./f)
    return (phases, times, trials)

def createJobDir(options):
    # Create job directory in options.output_dir/jobXXXX
    outputDir = options.output_dir + '/job{0:04}'.format(options.job_num) +'/'
    ec = system('mkdir ' + outputDir)
    if ec != 0:
        print "Could not create output directory: " + outputDir
        ec = system('ls ' + outputDir)
        if ec == 0:
            print "But it can be listed --> continuing!"
        else:
            print "And it cannot be listed. Check your permissions and rerun!"
            exit(1)
    return outputDir

def phaseCWT(sig, Tph, dt, maxF, dF=2):
    '''
    Calculate Morlet wavelet transform of a signal, but as a function of
    phase. Unaligned phase at the end will be discarded, and ph(t=0) must be 0,
    i.e. no phase shifts!
    '''
    n_ph = Tph/dt
    N = len(sig)
    q_ph = np.floor(N/n_ph)

    minF = 1./(len(sig)/2 * Morlet.fourierwl * dt)
    F = np.linspace(minF, maxF, (maxF-minF)/dF+1)
    scales = 1/F * 1/Morlet.fourierwl * 1/dt

    w = Morlet(sig, scales, scaling='direct')
    w_cwt_ph = np.ndarray((w.nscale, n_ph))

    #import pdb; pdb.set_trace()
    for sc_it in xrange(w.nscale):
        w_ph = np.reshape(np.abs(w.cwt[sc_it, :][0:q_ph*n_ph])**2, (q_ph, n_ph))
        w_cwt_ph[sc_it, :] = np.mean(w_ph, 0)

    sig_ph = np.reshape(sig[0:q_ph*n_ph], (q_ph, n_ph))
    phases = 1. * np.arange(n_ph) / n_ph * 2*np.pi - np.pi
    return phases, w_cwt_ph, 1./(w.scales*w.fourierwl*dt), sig_ph



def CWT(sig, dt, maxF, dF=2):
    '''
    Calculate a Morlet wavelet transfrom of a signal.
    '''
    N = len(sig)

    minF = 1./(len(sig)/2 * Morlet.fourierwl * dt)
    F = np.linspace(minF, maxF, (maxF-minF)/dF+1)
    scales = 1/F * 1/Morlet.fourierwl * 1/dt

    w = Morlet(sig, scales, scaling='direct')
    return np.abs(w.cwt)**2, 1./(w.scales*w.fourierwl*dt)


