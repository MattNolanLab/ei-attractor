#!/usr/bin/env python
#
#   plot_gaussian_err.py
#
#   Plot the error (distance from zero) in a random walk simulation.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import numpy as np

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

dt    = 0.001 # s
T     = 300
sigma = 3.16 # Misestimation of velocity cm/s/sqrt(dt)
NTrials = 100

figure(figsize=(12, 6))
subplot(2, 1, 1)
nt = T/dt + 1
v = np.random.randn(nt, NTrials)*np.sqrt(dt)
s = sigma*np.cumsum(v, axis=0)
times = np.arange(nt)*dt
plot(times, np.abs(s), '-')
ylabel('s (cm)')

subplot(2, 1, 2)
plot(times, np.std(s, axis=1), '-')
ylabel('$\sigma_{s}$ (cm)')
xlabel('Time (s)')
tight_layout()

savefig('error.png', dpi=150)

