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
import scipy
from scipy.signal import *
from numpy.fft.fftpack import fft
from Wavelets import Morlet

from matplotlib.pyplot import *

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



def splitSigToThetaCycles(sig, thetaT, dt):
    '''
    Take a 1D signal (np array) and rescale it to a 2D array, in which every row
    corresponds to one theta cycle.
    thetaT must be a multiple of dt
    The last, unaligned part of the signal will be discarded.

    Phase(sig, t=0) must be 0, no phase shifts!
    '''
    n_ph = thetaT / dt
    q_ph = len(sig) // n_ph
    return np.reshape(sig[0:q_ph*n_ph], (q_ph, n_ph))


def getChargeTheta(sig_theta_sliced, dt):
    '''
    For each theta cycle, find the total charge of synaptic current.
    Each row of sig_theta_sliced is one theta cycle
    '''
    return np.trapz(sig_theta_sliced, dx=dt, axis=1)



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


def createIgorSpikeRaster(spikes, yvals=None):
    '''
    spikes  row-wise 2d array of spike times
    '''
    if yvals is None:
        yvals = range(len(spikes)) + 1

    raster_x = np.ndarray((0))
    raster_y = np.ndarray((0))
    for it in range(len(spikes)):
        raster_x = np.hstack((raster_x, spikes[it]))
        raster_y = np.hstack((raster_y, np.zeros((len(spikes[it]))) + yvals[it]))

    return (raster_x, raster_y)

def fft_real_freq(sig, dt):
    '''
    Compute a DFT of a real signal and return an array of frequencies and
    Fourier coefficients
    '''
    S = fft(sig)
    S_F = np.linspace(0, 1, len(S)/2) / dt / 2.0

    return S_F, S[0:len(S_F)]


    
##############################################################################
#                         Image analysis functions
##############################################################################

def fitGaussian2D(sig, X, Y, C0, A0, mu0_x, mu0_y, sigma0):
    '''
    Fit a 2D Gaussian function to a 2D signal
    rateMap     A 2D array, population firing rate
    X, Y        Spatial locations, the same size as rateMap
    C0, A0, mu0_x, mu0_y, sigma0 
                Initial parameters for the Gaussian

        fun = C + A*exp(-||X - mu||^2 / (2*sigma^2))
    '''
    x0 = np.array([C0, A0, mu0_x, mu0_y, sigma0])
    fun = lambda x:  (x[0] + x[1] * np.exp( -((X - x[2])**2 + (Y - x[3])**2)/2./ x[4]**2 )) - sig
    #                  |      |                      |               |             |
    #                  C      A                    mu_x            mu_y          sigma

    xest,ierr = scipy.optimize.leastsq(fun, x0)
    return xest


def fitGaussianBump2D(rateMap, X, Y):
    '''
    Fit a Gaussian to a rate map, using least squares method.
    '''
    C0      = 0.0
    A0      = 1.0
    mu0_x   = 1.0
    mu0_y   = 1.0
    sigma0  = 1.0
    
    return fitGaussian2D(rateMap.flatten(), X.flatten(), Y.flatten(), C0, A0,
            mu0_x, mu0_y, sigma0)
    


##############################################################################
#                                   Tests
##############################################################################

if __name__ == "__main__":
    sz = (100, 100)

    X, Y = np.meshgrid(np.arange(sz[1]), np.arange(sz[0]))
    X -= 0.5 * sz[1]
    Y -= 0.5 * sz[0]

    C       = 0.0
    A       = 10.0
    mu_x    = -3.0
    mu_y    = 11.5
    sigma   = 100.0

    rateMap = C + A*np.exp(- ((X - mu_x)**2 + (Y - mu_y)**2) / (2*sigma**2))
    param_est = fitGaussianBump2D(rateMap, X, Y)
    Cest, Aest, muxest, muyest, sigmaest = param_est

    print param_est
    
    



