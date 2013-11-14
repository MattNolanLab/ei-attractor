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
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from plotting.global_defs import globalAxesSettings, createColorbar

linewidth=1

dx = 0.001
x0 = -0.5
x1 = 0.5

def plotWeights(ax, d, exc_profile, inh_profile, inh_const):
    plt.hold('on')
    ax = plt.gca()
    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ep, = plt.plot(d, exc_profile, linewidth=linewidth, color='red', label="E")
    ip, = plt.plot(d, inh_profile, linewidth=linewidth, color='blue', label="I")
    icp, = plt.plot(d, [inh_const]*len(d), ':', color='blue')
    plt.xlabel('Distance')
    plt.ylabel('G/G$_\mathrm{max}$')
    ax.yaxis.set_ticks([0, 1])
    ax.xaxis.set_ticks([x0, 0, x1])
    leg1 = ['E$\\rightarrow$I', 'I$\\rightarrow$E']
    leg2 = ['I$\\rightarrow$E uniform\nrandom']
    l1 = ax.legend([ep, ip], leg1, loc=(0.02, 1.0), frameon=False, fontsize='x-small',
            ncol=1)
    l2 = ax.legend([icp], leg2, loc=(0.45, 1.03), frameon=False, fontsize='x-small')
    plt.setp(l1.get_title(), fontsize='x-small')
    plt.setp(l2.get_title(), fontsize='x-small')
    ax.add_artist(l1)
    ax.margins(0.02)

    #arrow_clr='grey'
    #arrowprops = dict(
    #    arrowstyle = "->",
    #    linewidth=.5,
    #    color = arrow_clr,
    #    connectionstyle = "angle,angleA=0,angleB=90,rad=10")
    #
    #rnd_x, rnd_y = 0., inh_const
    #plt.annotate('Random uniform',
    #            (rnd_x, rnd_y), xytext=(0.4, 1.2), textcoords='axes fraction',
    #            arrowprops=arrowprops, ha='left', va='center',  size='small',
    #            color=arrow_clr, zorder=-1)
    


if (__name__ == "__main__"):
    from matplotlib import rc
    rc('pdf', fonttype=42)
    rc('mathtext', default='regular')

    plt.rcParams['font.size'] = 11

    d = np.arange(x0, x1+dx, dx)
    y_dim = np.sqrt(3)/2.0
    ES_pAMPA_mu = y_dim/2.0
    ES_pAMPA_sigma = 0.5/6
    ES_pGABA_sigma = 0.5/6
    ES_pGABA_const = 0.1
    shift = 0.1

    figsize = (2.5, 1.8)
    left    = 0.2
    bottom  = 0.3
    top     = 0.8
    right   = 0.85

    # Excitatory surround
    ES_exc_profile         = np.exp(-(np.abs(d) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
    ES_exc_profile_shifted = np.exp(-(np.abs(d - shift) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
    ES_inh_profile         = (1-ES_pGABA_const)*np.exp(-d**2/2./ES_pGABA_sigma**2) + ES_pGABA_const

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotWeights(ax, d, ES_exc_profile, ES_inh_profile, ES_pGABA_const)
    fileBase = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    fileName = "{0}/{1}_E_surr.pdf".format("output", fileBase)
    plt.savefig(fileName, transparent=True)

    # Inhibitory surround
    IS_pAMPA_sigma = 0.5/6
    IS_pGABA_mu    = y_dim/2.0
    IS_pGABA_sigma =  0.5/6
    IS_pGABA_const =  0.1
    IS_exc_profile = np.exp(-d**2/2./IS_pAMPA_sigma**2)
    IS_inh_profile = (1-IS_pGABA_const)*np.exp(-(np.abs(d) - IS_pGABA_mu)**2/2/IS_pGABA_sigma**2) + IS_pGABA_const

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotWeights(ax, d, IS_exc_profile, IS_inh_profile, IS_pGABA_const)
    fileBase = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    fileName = "{0}/{1}_I_surr.pdf".format("output", fileBase)
    plt.savefig(fileName, transparent=True)

    
