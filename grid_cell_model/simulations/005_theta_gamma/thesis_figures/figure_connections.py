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
import matplotlib.ticker as ti
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

pcolorFigsize     = (1.5, 1.5)
pcolorLeft        = 0.3
pcolorBottom      = 0.3
pcolorRight       = 0.92
pcolorTop         = 0.85
pcolorTransparent = True
##############################################################################


def plotHistogram(sp, neuronIdx, type, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', (1.75, 1.75))
    left        = kw.pop('left', 0.37)
    bottom      = kw.pop('bottom', 0.32)
    right       = kw.pop('right', 0.92)
    top         = kw.pop('top', 0.85)
    transparent = kw.pop('transparent', True)
    xlim        = kw.pop('xlim', None)
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

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.savefig(fname, dpi=150, transparent=transparent)



def plotOutgoing(sp, type, neuronIdx, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', pcolorFigsize)
    left        = kw.pop('left', pcolorLeft)
    bottom      = kw.pop('bottom', pcolorBottom)
    right       = kw.pop('right', pcolorRight)
    top         = kw.pop('top', pcolorTop)
    transparent = kw.pop('transparent', pcolorTransparent)

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
    pconn.plot2DWeightMatrix(conns, ax=ax, labelpad=1, **kw)
    ax.tick_params(which='both', length=0, pad=0)
    plt.savefig(fname, dpi=150, transparent=transparent)
    plt.close()


def plotIncoming(sp, type, neuronIdx, fname, trialNum=0, **kw):
    # kw arguments
    figsize     = kw.pop('figsize', pcolorFigsize)
    left        = kw.pop('left', pcolorLeft)
    bottom      = kw.pop('bottom', pcolorBottom)
    right       = kw.pop('right', pcolorRight)
    top         = kw.pop('top', pcolorTop)
    transparent = kw.pop('transparent', pcolorTransparent)

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
    pconn.plot2DWeightMatrix(conns, ax=ax, **kw)
    ax.tick_params(which='both', length=0, pad=0)
    plt.savefig(fname, dpi=150, transparent=transparent)
    plt.close()



def setFName(fname, outputRoot='output'):
    return "{0}/{1}".format(outputRoot, fname)

##############################################################################
EneuronIdx = 1937
IneuronIdx = 492

E_surrSp = JobTrialSpace2D(shape, E_surrRoot)
I_surrSp = JobTrialSpace2D(shape, I_surrRoot)

if (hists):
    # E-surround
    fname = setFName("conn_histogram_E_surr_E2I.pdf")
    plotHistogram(E_surrSp, IneuronIdx, "E2I", fname,
            locators=dict(
                x_major=ti.MultipleLocator(0.4),
                x_minor=ti.MultipleLocator(0.1)))

    fname = setFName("conn_histogram_E_surr_I2E.pdf")
    plotHistogram(E_surrSp, EneuronIdx, "I2E", fname,
            ylabel='', yticks=False,
            locators=dict(
                x_major=ti.MultipleLocator(1),
                x_minor=ti.MultipleLocator(0.5)))

    # I-surround
    fname = setFName("conn_histogram_I_surr_E2I.pdf")
    plotHistogram(I_surrSp, IneuronIdx, "E2I", fname,
            locators=dict(
                x_major=ti.MultipleLocator(1),
                x_minor=ti.MultipleLocator(0.25)))

    fname = setFName("conn_histogram_I_surr_I2E.pdf")
    plotHistogram(I_surrSp, EneuronIdx, "I2E", fname,
            ylabel='', yticks=False,
            locators=dict(
                x_major=ti.MultipleLocator(0.5),
                x_minor=ti.MultipleLocator(0.25)))




if (weights):
    # E-surround ######
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fname = setFName("conn_pcolor_E_surr_out_E2I.pdf")
    plotOutgoing(E_surrSp, "E", EneuronIdx, fname,
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_E_surr_out_I2E.pdf")
    plotOutgoing(E_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fname = setFName("conn_pcolor_E_surr_in_E2I.pdf")
    plotIncoming(E_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_E_surr_in_I2E.pdf")
    plotIncoming(E_surrSp, "E", EneuronIdx, fname,
            ylabel='',
            xlabel='')


    # I-surround ######
    fname = setFName("conn_pcolor_I_surr_out_E2I.pdf")
    plotOutgoing(I_surrSp, "E", EneuronIdx, fname,
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_I_surr_out_I2E.pdf")
    plotOutgoing(I_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fname = setFName("conn_pcolor_I_surr_in_E2I.pdf")
    plotIncoming(I_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # I-->E
    fname = setFName("conn_pcolor_I_surr_in_I2E.pdf")
    plotIncoming(I_surrSp, "E", EneuronIdx, fname,
            ylabel='',
            xlabel='')



