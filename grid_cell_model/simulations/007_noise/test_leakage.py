#
#   test_leakage.py
#
#   Test the "spectral leakage" hypothesis on peaks which are multiple of theta
#   frequency
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
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, xlabel, ylabel, show, xlim, figure
from matplotlib.mlab import psd, window_hanning, window_none
from scipy.signal    import correlate


F = [8., 48. , 56.  , 64.]
A = [0.,   .25,   .57,  0.26]
#P = [0., np.pi/4., np.pi/2., np.pi/4.*3.]
noiseA = 0.0

dt = 0.1e-3
T  = 1
t  = np.arange(0., T, dt)

sig = np.ndarray((len(t), ))
sig[:] = .0
for i in range(len(F)):
    sig += A[i] * np.cos(2 * np.pi * F[i] * t) # + P[i]) 
sig += noiseA * np.random.randn(len(sig))

figure(figsize=(10, 6))
plot(t, sig)
xlabel('Time (s)')
ylabel('Amplitude')


figure()
NFFT = len(sig)
#noverlap = int(NFFT*0.75)
noverlap = 0
win = window_none
Pxx, F = psd(sig - np.mean(sig), NFFT, Fs=1./dt, noverlap=noverlap,
        window=win)
plot(F, Pxx)
xlabel('Frequency (Hz)')
ylabel('Power (1/Hz)')
xlim([0, 100])

figure()
plt.plot(t, correlate(sig, sig)[len(sig) - 1:])

show()
