#!/usr/bin/env python
#
#   fig_conn_func.py
#
#   Print synaptic connection profiles.
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
import sys
import os

from matplotlib.pyplot  import figure, subplot, hold, plot, xlabel, ylabel, \
        gca, legend, setp, xlim, annotate, tight_layout, rcParams, savefig

from plotting.global_defs import globalAxesSettings, createColorbar

linewidth=1

dx = 0.001
x0 = -0.5
x1 = 0.5
d = np.arange(x0, x1+dx, dx)

y_dim = np.sqrt(3)/2.0
pAMPA_mu = y_dim/2.0
pAMPA_sigma = 0.5/6
pGABA_mu    = y_dim/2.0
pGABA_sigma = 0.5/6
pGABA_const = 0.1

shift = 0.1

# Excitatory surround
ES_exc_profile         = np.exp(-(np.abs(d) - pAMPA_mu)**2/2/pAMPA_sigma**2)
ES_exc_profile_shifted = np.exp(-(np.abs(d - shift) - pAMPA_mu)**2/2/pAMPA_sigma**2)
ES_inh_profile         = (1-pGABA_const)*np.exp(-d**2/2./pGABA_sigma**2) + pGABA_const
# Inhibitory surround
IS_exc_profile         = np.exp(-d**2/2./pAMPA_sigma**2)
IS_inh_profile         = (1-pGABA_const)*np.exp(-(d - pGABA_mu)**2/2/pGABA_sigma**2) + pGABA_const


def plotWeights(ax):
    hold('on')
    ax = gca()
    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plot(d, ES_exc_profile, linewidth=linewidth, color='red', label="E")
    plot(d, [pGABA_const]*len(d), ':', color='blue')
    plot(d, ES_inh_profile, linewidth=linewidth, color='blue', label="I")
    xlabel('Distance')
    ylabel('G/G$_\mathrm{max}$')
    ax.yaxis.set_ticks([0, 1])
    ax.xaxis.set_ticks([x0, 0, x1])
    #legend(bbox_to_anchor=(0., 1.05, 1., 1.05), ncol=2, loc=3, mode='expand', borderaxespad=0.)
    #setp(ax.get_legend().get_texts(), fontsize='small')
    #xlim([x0, x1])
    ax.margins(0.02)

    arrow_clr='grey'
    arrowprops = dict(
        arrowstyle = "->",
        linewidth=.5,
        color = arrow_clr,
        connectionstyle = "angle,angleA=0,angleB=90,rad=10")
    
    rnd_x, rnd_y = 0., pGABA_const
    annotate('Random uniform',
                (rnd_x, rnd_y), xytext=(0.4, 1.2), textcoords='axes fraction',
                arrowprops=arrowprops, ha='left', va='center',  size='small',
                color=arrow_clr, zorder=-1)
    


if (__name__ == "__main__"):
    from matplotlib import rc
    rc('pdf', fonttype=42)
    rc('mathtext', default='regular')

    rcParams['font.size'] = 11
    figSize = (2.5, 1.8)
    fig = figure(figsize=figSize)
    ax = subplot(111)
    fig.subplots_adjust(left=0.2, top=0.8, bottom=0.3, right=0.85)
    plotWeights(ax)
    #tight_layout(rect=(0., 0., 0.99, 0.9))
    
    
    if (len(sys.argv) > 1):
        fileName = sys.argv[1]
    else:
        fileBase = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        fileName = "{0}/{1}.pdf".format("output", fileBase)
    savefig(fileName, transparent=True)
    
