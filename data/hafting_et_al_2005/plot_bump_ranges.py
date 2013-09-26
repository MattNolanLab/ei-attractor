#!/usr/bin/env python
#
#   plot_bump_ranges.py
#
#   Plot velocity ranges and bump ranges needed to cover a specified grid field
#   spacing.
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
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, \
        ScalarFormatter, MaxNLocator

from plotting.grids import gridScaleBar
from plotting.global_defs import globalAxesSettings

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

velFileName = 'rat_trajectory_lowpass.mat'
data = loadmat(velFileName)

dt = float(data['dt'])
pos_x = data['pos_x'].flatten()
pos_y = data['pos_y'].flatten()


# Positions 
figure()
plot(pos_x, pos_y)
hold('on')
axis('off')
axis('scaled')
title('Positions: low pass filtered', va='bottom')
gridScaleBar(50, True, gca())
savefig('positions_lowpass.pdf')


def setAxes(ax):
    globalAxesSettings(ax)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    f = ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits((0, 3))
    ax.yaxis.set_major_formatter(f)

# histograms of velocities
hist_bins = 40
grid_lambda_all = [40, 50, 60.0] # grid spacing: cm
Nx = 34.0          # size of neural sheet X (neurons)
vel_x = np.diff(pos_x)/dt
vel_y = np.diff(pos_y)/dt
vel = np.hstack((vel_x, vel_y))
figure(figsize=(6.0, 3.0))
subplot(1, 2, 1)
ax = gca()
setAxes(ax)
hist(np.abs(vel), bins=hist_bins, normed=True, histtype='bar')
xlabel('Animal speed (cm/s)')
ylabel('p($\cdot$)')
title('A', fontweight='bold', x=-0.35, y=1.1, ha='left')
subplot(1, 2, 2)
ax = gca()
setAxes(ax)
leg = []
for grid_lambda in grid_lambda_all:
    bump_s = Nx / grid_lambda * np.abs(vel)
    hist(bump_s, bins=hist_bins, normed=True, histtype='step')
    leg.append('{0}'.format(int(grid_lambda)))
l = legend(leg, title='$\lambda_{grid}$ (cm)', fontsize='small', frameon=False)
plt.setp(l.get_title(), fontsize='small')
xlabel('Bump speed (neurons/s)')
title('B', fontweight='bold', x=-0.2, y=1.1, ha='left')
tight_layout(w_pad=1, pad=0.2)
savefig('velocity_histograms.pdf')
savefig('velocity_histograms.eps')
