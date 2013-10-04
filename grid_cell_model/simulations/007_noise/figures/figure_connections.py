#!/usr/bin/env python
#
#   figure_connections.py
#
#   Noise publication figures: connection weight figures.
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
from matplotlib.ticker     import MultipleLocator, AutoMinorLocator, \
        MaxNLocator
from matplotlib.transforms import Bbox

from EI_plotting            import computeYX
from parameters.param_space import JobTrialSpace2D, DataSpace
from plotting.global_defs   import globalAxesSettings
from plotting.connections   import plot2DWeightMatrix

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11


DS = DataSpace

##############################################################################
connDataRoot= 'output_local/even_spacing/connections'
shape = (1, 31)
iterList  = ['g_AMPA_total', 'g_GABA_total']

hists   = 1
weights = 0
##############################################################################

def plotConnHistogram(val, **kw):
    # keyword arguments
    kw['bins']      = kw.get('bins', 20)
    #kw['edgecolor'] = kw.get('edgecolor', 'none')
    ax              = kw.get('ax', plt.gca())
    xlabel          = kw.pop('xlabel', 'g (nS)')
    ylabel          = kw.pop('ylabel', 'Count')
    title           = kw.pop('title', '')
    locators        = kw.pop('locator', {})

    globalAxesSettings(ax)
    ax.hist(val, **kw)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # tick formatting
    x_major = locators.get('x_major', MultipleLocator(1))
    x_minor = locators.get('x_minor', AutoMinorLocator(2))
    y_major = locators.get('y_major', MaxNLocator(4))
    y_minor = locators.get('y_minor', AutoMinorLocator(2))
    ax.xaxis.set_major_locator(x_major)
    ax.yaxis.set_major_locator(y_major)
    ax.xaxis.set_minor_locator(x_minor)
    ax.yaxis.set_minor_locator(y_minor)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


def plotEToI(sp, gIdx, neuronIdx, trialNum=0, **kw):
    gE, gI = computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_IE']
    conns  = M[neuronIdx, :]
    ax = plotConnHistogram(conns,
            title='I cell', **kw)
    annG = gE[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_E$ = {0} (nS)'.format(annG)
    ax.text(0.95, 0.95, ann, ha='right', va='top', fontsize='small',
            transform=ax.transAxes)

def plotIToE(sp, gIdx, neuronIdx, trialNum=0, **kw):
    gE, gI = computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]
    ax = plotConnHistogram(conns,
            title='E cell', **kw)
    annG = gI[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_I$ = {0} (nS)'.format(annG)
    ax.text(0.95, 0.95, ann, ha='right', va='top', fontsize='small',
            transform=ax.transAxes)


def plotOutgoing(sp, gIdx, type, neuronIdx, trialNum=0, **kw):
    Nx     = None
    Ni     = None

    data = sp[0][gIdx][trialNum].data
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
    plot2DWeightMatrix(conns, **kw)

def plotIncoming(sp, gIdx, type, neuronIdx, trialNum=0, **kw):
    Nx = None
    Ni = None

    data = sp[0][gIdx][trialNum].data
    if (type == 'I'):
        var = 'g_IE'
        Nx = DS.getNetParam(data, 'Ne_x')
        Ny = DS.getNetParam(data, 'Ne_y')
        kw['title'] = 'E cells $\\rightarrow$ I cell'
        
    elif (type == 'E'):
        var = 'g_EI'
        Nx = DS.getNetParam(data, 'Ni_x')
        Ny = DS.getNetParam(data, 'Ni_y')
        kw['title'] = 'I cells $\\rightarrow$ E cell'

    conns = np.reshape(data[var][neuronIdx, :], (Ny, Nx))
    plot2DWeightMatrix(conns, **kw)





##############################################################################
gIdx = 15
neuronIdx = 0

figSize = (1.75, 1.75)
left   = 0.35
bottom = 0.32
right  = 0.95
top    = 0.85
transparent = True

sp = JobTrialSpace2D(shape, connDataRoot)

if (hists):
    fig = plt.figure('E2I', figsize=figSize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotEToI(sp, gIdx, neuronIdx)
    #fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=0)
    fname = "figure_connections_E2I.pdf"
    plt.savefig(fname, dpi=300, transparent=transparent)


    fig = plt.figure('I2E', figsize=figSize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotIToE(sp, gIdx, neuronIdx, ylabel='')
    #fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=0)
    fname = "figure_connections_I2E.pdf"
    plt.savefig(fname, dpi=300, transparent=transparent)


if (weights):
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fig = plt.figure('g_out_E2I', figsize=figSize)
    plotOutgoing(sp, gIdx, "E", neuronIdx)
    fig.tight_layout()
    fname = "figure_connections_pcolor_out_E2I.png"
    plt.savefig(fname, dpi=300, transparent=False)

    # I-->E
    fig = plt.figure('g_out_I2E', figsize=figSize)
    plotOutgoing(sp, gIdx, "I", neuronIdx)
    fig.tight_layout()
    fname = "figure_connections_pcolor_out_I2E.png"
    plt.savefig(fname, dpi=300, transparent=False)


    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fig = plt.figure('g_in_E2I', figsize=figSize)
    plotIncoming(sp, gIdx, "I", neuronIdx)
    fig.tight_layout()
    fname = "figure_connections_pcolor_in_E2I.png"
    plt.savefig(fname, dpi=300, transparent=False)

    # I-->E
    fig = plt.figure('g_in_I2E', figsize=figSize)
    plotIncoming(sp, gIdx, "E", neuronIdx)
    fig.tight_layout()
    fname = "figure_connections_pcolor_in_I2E.png"
    plt.savefig(fname, dpi=300, transparent=False)



