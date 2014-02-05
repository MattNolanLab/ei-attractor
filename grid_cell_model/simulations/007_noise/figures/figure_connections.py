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

import plotting.connections as pconn
from EI_plotting            import aggregate as aggr
from parameters.param_space import JobTrialSpace2D, DataSpace
from plotting.global_defs   import globalAxesSettings
from submitting import flagparse

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11
outputDir = "panels"


DS = DataSpace

##############################################################################

connDataRoot= 'output_local/even_spacing/connections'
shape = (1, 31)
iterList  = ['g_AMPA_total', 'g_GABA_total']

parser = flagparse.FlagParser()
parser.add_flag('--hists')
parser.add_flag('-w', '--weights')
args = parser.parse_args()

##############################################################################

def plotEToI(sp, gIdx, neuronIdx, trialNum=0, **kw):
    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_IE']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
            title='I cell', **kw)
    annG = gE[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_E$ = {0} (nS)'.format(annG)
    ax.text(0.95, 0.95, ann, ha='right', va='top', fontsize='small',
            transform=ax.transAxes)

def plotIToE(sp, gIdx, neuronIdx, trialNum=0, **kw):
    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
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
    pconn.plot2DWeightMatrix(conns, **kw)

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
    pconn.plot2DWeightMatrix(conns, **kw)





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

if args.hists or args.all:
    fig = plt.figure('E2I', figsize=figSize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotEToI(sp, gIdx, neuronIdx)
    #fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=0)
    fname = outputDir + "/figure_connections_E2I.pdf"
    plt.savefig(fname, dpi=300, transparent=transparent)


    fig = plt.figure('I2E', figsize=figSize)
    ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
    plotIToE(sp, gIdx, neuronIdx, ylabel='')
    #fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=0)
    fname = outputDir + "/figure_connections_I2E.pdf"
    plt.savefig(fname, dpi=300, transparent=transparent)


if args.weights or args.all:
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fig = plt.figure('g_out_E2I', figsize=figSize)
    plotOutgoing(sp, gIdx, "E", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_out_E2I.png"
    plt.savefig(fname, dpi=300, transparent=False)

    # I-->E
    fig = plt.figure('g_out_I2E', figsize=figSize)
    plotOutgoing(sp, gIdx, "I", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_out_I2E.png"
    plt.savefig(fname, dpi=300, transparent=False)


    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fig = plt.figure('g_in_E2I', figsize=figSize)
    plotIncoming(sp, gIdx, "I", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_in_E2I.png"
    plt.savefig(fname, dpi=300, transparent=False)

    # I-->E
    fig = plt.figure('g_in_I2E', figsize=figSize)
    plotIncoming(sp, gIdx, "E", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_in_I2E.png"
    plt.savefig(fname, dpi=300, transparent=False)



