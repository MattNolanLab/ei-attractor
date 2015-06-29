'''Signal analysis tools.

The :mod:`analysis.signal` module contains functions and classes for continuous
signal analysis. This can be e.g. filtering, slicing, correlation analysis,
up/down-sampling, etc.

.. note::
    Some of these functions are re-implemented and with much better
    documentation in the gridcells_ package.

.. currentmodule:: grid_cell_model.analysis.signal

.. _gridcells: http://gridcells.readthedocs.org

Functions
---------

.. autosummary::

    butterHighPass
    butterBandPass
    spikePhaseTrialRaster
    splitSigToThetaCycles
    getChargeTheta
    phaseCWT
    CWT
    fft_real_freq
    relativePower
    maxPowerFrequency
    localExtrema
    globalExtremum
    relativePeakHeight
    downSample
    sliceSignal
'''

import numpy as np
import scipy.signal

from numpy.fft.fftpack import fft
from Wavelets import Morlet

__all__ = ['butterHighPass', 'butterBandPass', 'spikePhaseTrialRaster',
           'splitSigToThetaCycles', 'getChargeTheta', 'phaseCWT', 'CWT',
           'fft_real_freq', 'relativePower', 'maxPowerFrequency',
           'localExtrema', 'globalExtremum', 'relativePeakHeight',
           'downSample', 'sliceSignal']


def butterHighPass(sig, dt, f_pass):
    nyq_f = 1./dt/2
    norm_f_pass = f_pass/nyq_f

    # Low pass filter
    b, a = scipy.signal.butter(3, norm_f_pass, btype='high')
    return scipy.signal.filtfilt(b, a, sig)


def butterBandPass(sig, dt, f_start, f_stop):
    '''Band pass filter a signal, with f_start and f_stop frequencies'''
    nyq_f = 1./dt/2
    norm_f_start = f_start/ nyq_f
    norm_f_stop  = f_stop / nyq_f
    b, a = scipy.signal.butter(3, [norm_f_start, norm_f_stop], btype='band')
    return scipy.signal.filtfilt(b, a, sig)


def spikePhaseTrialRaster(spikeTimes, f, start_t=0):
    '''Here assuming that phase(t=0) = 0'''
    spikeTimes -= start_t
    trials = np.floor(f*spikeTimes)
    phases = np.mod(2*np.pi*f*spikeTimes, 2*np.pi)
    times  = np.mod(spikeTimes, 1./f)
    return (phases, times, trials)




## Take a 1D signal and rescale it to signals of individual theta cycles.
#
# Each row of the result contains one theta cycle and it is assumed that theta
# is generated continuously.
#
# The last, unaligned part of the signal will be discarded.
#
# Phase(sig, t=0) must be 0, no phase shifts!
#
# @param sig   Signal with dt.
# @param thetaT Theta period. MUST be a multiple of dt
# @param dt     Time resolution of the signal.
def splitSigToThetaCycles(sig, thetaT, dt):
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



## Compute a DFT of a real signal and return an array of frequencies and
# Fourier coefficients.
#
# @param sig    The signal
# @param dt     Sampling rate of the signal
# @return A tuple (F, sig~). F is an array of frequencies, sig~ is the fourier
#         transform
#
def fft_real_freq(sig, dt):
    S = fft(sig)
    S_F = np.linspace(0, 1, len(S)/2) / dt / 2.0

    return S_F, S[0:len(S_F)]





## Compute power from FFT data in a specified frequency range, relative to the
# total power.
#
# This function will throw an error if the desired frequency range is out of the
# range of the actual signal.
#
# @param Pxx     A Power spectral density vector
# @param F       Frequencies corresponding to Pxx (Hz).
# @param Frange  A tuple containing the frequency range (Hz).
# @return Relative power in the specified freqency range.
#
def relativePower(Pxx, F, Frange):
    Fidx = np.nonzero(np.logical_and(F >= Frange[0], F <= Frange[1]))[0]
    rangeP = sum(Pxx[Fidx])
    return rangeP / sum(Pxx)




# Get frequency with maximum power
#
# @param Pxx     Power spectral density of the signal.
# @param F       A corresponding array of frequencies
# @param Frange  A tuple containing frequency range to restrict the analysis to.
#
# @return An index to F, the frequency with maximum power.
#
def maxPowerFrequency(Pxx, F, Frange=None):
    if (Frange == None):
        return F[np.argmax(Pxx)]
    else:
        range = np.logical_and(F >= Frange[0], F <= Frange[1])
        return F[range][np.argmax(Pxx[range])]



###############################################################################
#                           Extrema analysis
###############################################################################
def localExtrema(sig):
    '''
    Find all local extrema using the derivative approach.

    Parameters
    ----------
    sig : numpy.ndarray
        A 1D numpy array

    output : (numpy.ndarray, numpy.ndarray)
        A pair (idx, types) containing the positions of local extrema iniside
        ``sig`` and the type of the extrema:

            * type > 0 means local maximum
            * type < 0 is local minimum
    '''
    sz = len(sig)
    szDiff = sz - 1
    der = np.diff(sig)
    der0 = (der[0:szDiff - 1] * der[1:szDiff]) < 0.
    ext_idx = np.nonzero(der0)[0]
    dder = np.diff(der)[ext_idx]
    ext_idx += 1    # Correction for a peak position
    ext_t = np.ndarray((dder.size, ), dtype=int)
    ext_t[dder < 0] = 1
    ext_t[dder > 0] = -1
    return (ext_idx, ext_t)



def globalExtremum(sig, func):
    '''
    Return global maximum of a signal.

    ``sig``
        A 1D array, in which case a global extremum will be returned as a single
        value, or a 2D array, for which an array of extrema will be returned,
        one extremum for **each row**. Numpy arrays accepted only.

    ``func``
        Numpy extremum function. Must accept the 'axis' parameter
    '''
    shpLen = len(sig.shape)
    if (shpLen == 1):
        return func(sig)
    elif (shpLen == 2):
        return func(sig, axis=1)
    else:
        raise TypeError("signal must be either 1D or a 2D numpy array!")


def relativePeakHeight(localExtrema, cmpFun):
    sz = len(localExtrema)
    if (sz == 0):
        raise TypeError("Cannot compute relative peak heights in an empty array.")
    elif (sz == 1):
        raise TypeError("Cannot compute relative peak heights in an array with only one element")

    res = np.ndarray((sz, ))
    cmp = np.ndarray((2, sz-2))

    hr = np.abs(localExtrema[0:sz-1] - localExtrema[1:])
    hl = np.abs(localExtrema[1:]     - localExtrema[0:sz-1])


    cmp[0, :] = hr[1:]
    cmp[1, :] = hl[0:hl.size - 1]

    res[0]  = hr[0]
    res[-1] = hl[-1]
    res[1:sz-1] = cmpFun(cmp, axis=0)

    #import pdb; pdb.set_trace()

    return res



###############################################################################
#                         (Re-)Sampling signals
###############################################################################

def downSample(sig, factor, X=None):
    '''
    Downsample a signal by a specified factor.

    Parameters
    ----------
    sig : numpy.ndarray
        Signal to down-sample.

    factor : int
        Down-sampling factor. If the factor is not int, it will be converted to
        int.

    X : numpy.ndarray, optional
        Optional X value of the signal (i.e. time, space) that will be
        downsampled as well. If ``None``, [0, len(sig)-1] will be used.

    output : numpy.ndarray
        A tuple (sig_d, X_d), where
        sig_d is the downsampled signal and X_d is the down-sampled X
        coordinate of the signal. The size of the output will depend on whether
        the size of ``sig`` is an integer multiple of ``factor``.
    '''
    if (X is None):
        X = np.arange(len(sig))
    idx = np.arange(0, len(sig), 10)
    return (sig[idx], X[idx])


###############################################################################
#                                     Other
###############################################################################
def sliceSignal(t, sig, tStart, tEnd):
    idx = np.logical_and(t >= tStart, t <= tEnd)
    return t[idx], sig[idx], idx

