#
#   bumps.py
#
#   Functions for plotting attractor bump related data.
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
import matplotlib.pyplot as mpl
import numpy as np
from matplotlib.pyplot  import plot, xlabel, ylabel, legend, ylim, \
        tight_layout, axis, title, pcolormesh, colorbar, hold, subplot
from matplotlib.ticker  import MaxNLocator, LinearLocator
from analysis.spikes    import torusPopulationVector

from global_defs        import globalAxesSettings, createColorbar


def bumpPosition(spikes, sheetSize, tstart, tend, dt, winLen, units="s"):
    '''
    Print a plot of bump position
    '''
    Ne_x = sheetSize[0]
    Ne_y = sheetSize[1]
    (pos, bumpPos_times) = torusPopulationVector(
            spikes,
            sheetSize,
            tstart,
            tend,
            dt,
            winLen)
    plot(bumpPos_times, pos)
    xlabel('Time (%s)' % units)
    ylabel('Bump position (neurons)')
    legend(['X', 'Y'])
    ylim([-Ne_x/2 - 5, Ne_x/2 + 5])
    tight_layout()



def torusFiringRate(rateMap, labelx, labely=None, titleStr=""):
    '''Firing rate of cells on the twisted torus'''
    if (labely is None):
        labely = labelx
    pcolormesh(rateMap)
    globalAxesSettings(mpl.gca())
    locx = LinearLocator(2)
    mpl.gca().xaxis.set_ticks([0, len(rateMap[0])])
    mpl.gca().yaxis.set_ticks([0, len(rateMap)])
    xlabel(labelx)
    ylabel(labely)
    createColorbar(mpl.gca(), data=rateMap, label='Firing rate (Hz)')
    axis('tight')
    title(titleStr, va='bottom')
    tight_layout()


def flatFiringRate(FR, times, labely=None, labelx=None, units="ms",
        titleStr="", nticks=2):
    globalAxesSettings(mpl.gca())
    N = len(FR) # no. of rows == no. of neurons
    T, N_id = np.meshgrid(times, np.arange(N))
    pcolormesh(T, N_id,  FR)
    loc = LinearLocator(nticks)
    loc.tick_values(0, N-1)
    mpl.gca().yaxis.set_major_locator(loc)
    if (labelx is None):
        xlabel("Time (%s)" % units)
    elif (labelx != ""):
        xlabel(labelx)

    if (labely is None):
        ylabel("Neuron #")
    elif (labely != ""):
        ylabel(labely)

    if (titleStr != ""):
        title(titleStr)

    axis('tight')
    createColorbar(mpl.gca(), data=FR, label='Firing rate (Hz)')
    tight_layout()


