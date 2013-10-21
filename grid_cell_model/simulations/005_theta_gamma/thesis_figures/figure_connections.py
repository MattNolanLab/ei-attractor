#!/usr/bin/env python
#
#   figure_connections.py
#
#   PhD. thesis figures: connection weight figures.
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
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

import plotting.connections as pconn
from parameters.param_space import JobTrialSpace2D, DataSpace

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

DS = DataSpace

##############################################################################
E_surrRoot = 'data/connections/E_surround'
I_surrRoot = 'data/connections/I_surround'
outputRoor = 'output'
shape = (1, 1)

hists   = 1
weights = 1
##############################################################################


def plotHistogram(sp, neuronIdx, type, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', (1.75, 1.75))
    left        = kw.pop('left', 0.37)
    bottom      = kw.pop('bottom', 0.32)
    right       = kw.pop('right', 0.92)
    top         = kw.pop('top', 0.85)
    transparent = kw.pop('transparent', True)
    ylim        = kw.pop('ylim', [0, 900])
    xticks      = kw.pop('xticks', True)
    yticks      = kw.pop('yticks', True)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))

    if (type == 'E2I'):
        dataVar = 'g_IE'
        title   = "I cell"
    elif (type == 'I2E'):
        dataVar = 'g_EI'
        title   = "E cell"
    else:
        raise ValueError()

    M      = sp[0][0][trialNum].data[dataVar]
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns, title=title, ax=ax, **kw)

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    ax.set_ylim(ylim)
    plt.savefig(fname, dpi=150, transparent=transparent)



def plotOutgoing(sp, type, neuronIdx, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', (1.75, 1.6))
    left        = kw.pop('left', 0.35)
    bottom      = kw.pop('bottom', 0.3)
    right       = kw.pop('right', 0.9)
    top         = kw.pop('top', 0.85)
    transparent = kw.pop('transparent', True)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))

    Nx = None
    Ni = None

    data = sp[0][0][trialNum].data
    if (type == 'E'):
        var = 'g_IE'
        Nx = DS.getNetParam(data, 'Ni_x')
        Ny = DS.getNetParam(data, 'Ni_y')
        kw['title'] = 'E cell $\\rightarrow$ I cells'
        
    elif (type == 'I'):
        var = 'g_EI'
        Nx = DS.getNetParam(data, 'Ne_x')
        Ny = DS.getNetParam(data, 'Ne_y')
        kw['title'] = 'I cell $\\rightarrow$ E cells'

    conns = np.reshape(data[var][:, neuronIdx], (Ny, Nx))
    pconn.plot2DWeightMatrix(conns, **kw)
    plt.savefig(fname, dpi=150, transparent=transparent)
    plt.close()


def plotIncoming(sp, type, neuronIdx, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', (1.75, 1.6))
    left        = kw.pop('left', 0.35)
    bottom      = kw.pop('bottom', 0.3)
    right       = kw.pop('right', 0.9)
    top         = kw.pop('top', 0.85)
    transparent = kw.pop('transparent', True)

    fig = plt.figure(figsize=figsize)
    ax = plt.gcf().add_axes(Bbox.from_extents(left, bottom, right, top))

    Nx = None
    Ni = None

    data = sp[0][0][trialNum].data
    if (type == 'I'):
        var = 'g_IE'
        Nx = DS.getNetParam(data, 'Ne_x')
        Ny = DS.getNetParam(data, 'Ne_y')
        kw['title'] = 'E cells$\\rightarrow$I cell'
        
    elif (type == 'E'):
        var = 'g_EI'
        Nx = DS.getNetParam(data, 'Ni_x')
        Ny = DS.getNetParam(data, 'Ni_y')
        kw['title'] = 'I cells$\\rightarrow$E cell'

    conns = np.reshape(data[var][neuronIdx, :], (Ny, Nx))
    pconn.plot2DWeightMatrix(conns, **kw)
    plt.savefig(fname, dpi=150, transparent=transparent)
    plt.close()



def setFName(fname, outputRoot='output'):
    return "{0}/{1}".format(outputRoot, fname)

##############################################################################
neuronIdx = 0

E_surrSp = JobTrialSpace2D(shape, E_surrRoot)

if (hists):
    fname = setFName("conn_histogram_E_surr_E2I.pdf")
    plotHistogram(E_surrSp, neuronIdx, "E2I", fname)

    fname = setFName("conn_histogram_E_surr_I2E.pdf")
    plotHistogram(E_surrSp, neuronIdx, "I2E", fname,
            ylabel='', yticks=False)




if (weights):
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fname = setFName("conn_pcolor_E_surr_out_E2I.pdf")
    plotOutgoing(E_surrSp, "E", neuronIdx, fname,
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_E_surr_out_I2E.pdf")
    plotOutgoing(E_surrSp, "I", neuronIdx, fname,
            ylabel='',
            xlabel='')

    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fname = setFName("conn_pcolor_E_surr_in_E2I.pdf")
    plotIncoming(E_surrSp, "I", neuronIdx, fname,
            ylabel='',
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_E_surr_in_I2E.pdf")
    plotIncoming(E_surrSp, "E", neuronIdx, fname,
            ylabel='',
            xlabel='')



