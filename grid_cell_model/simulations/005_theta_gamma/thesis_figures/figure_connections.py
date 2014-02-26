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

import settings as se
import plotting.connections as pconn
from parameters.param_space import JobTrialSpace2D, DataSpace

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

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
    ax          = kw.pop('ax', plt.gca())
    transparent = kw.pop('transparent', True)
    plotTitle   = kw.pop('plotTitle', True)
    xlim        = kw.pop('xlim', None)
    ylim        = kw.pop('ylim', [0, 900])
    xticks      = kw.pop('xticks', True)
    yticks      = kw.pop('yticks', True)

    if (type == 'E2I'):
        dataVar = 'g_IE'
        title   = "I cell"
    elif (type == 'I2E'):
        dataVar = 'g_EI'
        title   = "E cell"
    else:
        raise ValueError()

    if not plotTitle:
        title = ''

    M      = sp[0][0][trialNum].data[dataVar]
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns, title=title, ax=ax, **kw)

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    return ax



def plotBrokenHistogram(sp, neuronIdx, type, fname, trialNum=0,
        axBoundaries=[0, 0, 1, 1], axesProportions=(0.5, 0.5),
        bottomLimits=None, topLimits=None, **kw):
    # kw arguments
    plotTitle   = kw.pop('plotTitle', True)
    fig         = kw.pop('fig', plt.gcf())
    transparent = kw.pop('transparent', True)
    xlim        = kw.pop('xlim', None)

    left, bottom, right, top = axBoundaries
    h = top - bottom
    w = right - left
    hBottom = h*axesProportions[0]
    hTop = h*axesProportions[1]
    
    axBottom = fig.add_axes(Bbox.from_extents(left, bottom, right, bottom +
        hBottom))
    axTop = fig.add_axes(Bbox.from_extents(left, top - hTop, right, top),
            sharex=axBottom)

    if (type == 'E2I'):
        dataVar = 'g_IE'
        title   = "I cell"
    elif (type == 'I2E'):
        dataVar = 'g_EI'
        title   = "E cell"
    else:
        raise ValueError()

    if not plotTitle:
        title = ''

    M      = sp[0][0][trialNum].data[dataVar]
    conns  = M[neuronIdx, :]
    pconn.plotConnHistogram(conns, title=title, ax=axBottom, **kw)
    kw['ylabel'] = ''
    pconn.plotConnHistogram(conns, title=title, ax=axTop, **kw)

    #if (not xticks):
    #    ax.xaxis.set_ticklabels([])
    #if (not yticks):
    #    ax.yaxis.set_ticklabels([])

    axBottom.set_xlim(xlim)
    if xlim is not None:
        axBottom.set_xticks(xlim)
    axBottom.set_ylim(bottomLimits)
    axBottom.set_yticks(bottomLimits)
    axBottom.yaxis.set_minor_locator(ti.NullLocator())
    axTop.set_ylim(topLimits)
    axTop.set_yticks([topLimits[1]])
    axTop.xaxis.set_visible(False)
    axTop.spines['bottom'].set_visible(False)

    divLen = 0.07
    d = .015
    kwargs = dict(transform=fig.transFigure, color='k', clip_on=False)
    axBottom.plot((left-divLen*w, left+divLen*w), (bottom+hBottom + d,
        bottom+hBottom - d), **kwargs)
    axTop.plot((left-divLen*w, left+divLen*w), (top-hTop + d, top-hTop - d),
            **kwargs)

    return axBottom, axTop


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



##############################################################################
EneuronIdx = 1937
IneuronIdx = 492

E_surrSp = JobTrialSpace2D(shape, E_surrRoot)
I_surrSp = JobTrialSpace2D(shape, I_surrRoot)
histLeft   = 0.4
histBottom = 0.32
histRight  = 0.90
histTop    = 0.85
histFigsize = (1.8, 1.55)
axBoundaries = (histLeft, histBottom, histRight, histTop)

if (hists):
    # E-surround
    fname = se.setFName("conn_histogram_E_surr_E2I.pdf")
    fig = plt.figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotHistogram(E_surrSp, IneuronIdx, "E2I", fname,
            plotTitle=False,
            ax=ax,
            locators=dict(
                x_major=ti.MultipleLocator(0.4),
                x_minor=ti.MultipleLocator(0.1)),
            rwidth=0.8,
            linewidth=0)
    ax.set_xlabel(ax.xaxis.get_label_text(), labelpad=-5)
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

    fname = se.setFName("conn_histogram_E_surr_I2E.pdf")
    fig = plt.figure(figsize=histFigsize)
    axBottom, axTop = plotBrokenHistogram(E_surrSp, EneuronIdx, "I2E", fname,
            ylabel='', plotTitle=False,
            axBoundaries=axBoundaries,
            axesProportions=(0.75, 0.2),
            bottomLimits=(0, 60),
            topLimits=(800, 900),
            xlim=[0, 2.5],
            rwidth=0.8,
            linewidth=0)
    axBottom.set_xlabel(axBottom.xaxis.get_label_text(), labelpad=-5)
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()
            

    # I-surround
    fname = se.setFName("conn_histogram_I_surr_E2I.pdf")
    fig = plt.figure(figsize=histFigsize)
    axBottom, axTop = plotBrokenHistogram(I_surrSp, IneuronIdx, "E2I", fname,
            fig=fig,
            ylabel='',
            plotTitle=False,
            axBoundaries=axBoundaries,
            axesProportions=(0.75, 0.2),
            bottomLimits=(0, 100),
            topLimits=(800, 900),
            xlim=[0, 1.25],
            rwidth=0.8,
            linewidth=0)
    axBottom.set_xlabel(ax.xaxis.get_label_text(), labelpad=-5)
    fig.text(histLeft - 0.27, 0.5*(histBottom+histTop), 'Count', rotation=90,
            ha='center', va='center')
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

    fname = se.setFName("conn_histogram_I_surr_I2E.pdf")
    fig = plt.figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotHistogram(I_surrSp, IneuronIdx, "I2E", fname,
            ylabel='',
            plotTitle=False,
            ax=ax,
            xlim=[0, 0.8],
            ylim=[0, 150],
            locators=dict(
                x_major=ti.MultipleLocator(0.8),
                x_minor=ti.MultipleLocator(0.4)),
            rwidth=0.8,
            linewidth=0)
    ax.set_xlabel(ax.xaxis.get_label_text(), labelpad=-5)
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()




if (weights):
    # E-surround ######
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fname = se.setFName("conn_pcolor_E_surr_out_E2I.pdf")
    plotOutgoing(E_surrSp, "E", EneuronIdx, fname,
            xlabel='')

    # I-->E
    fname = se.setFName("conn_pcolor_E_surr_out_I2E.pdf")
    plotOutgoing(E_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fname = se.setFName("conn_pcolor_E_surr_in_E2I.pdf")
    plotIncoming(E_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # I-->E
    fname = se.setFName("conn_pcolor_E_surr_in_I2E.pdf")
    plotIncoming(E_surrSp, "E", EneuronIdx, fname,
            ylabel='',
            xlabel='')


    # I-surround ######
    fname = se.setFName("conn_pcolor_I_surr_out_E2I.pdf")
    plotOutgoing(I_surrSp, "E", EneuronIdx, fname,
            xlabel='')

    # I-->E
    fname = se.setFName("conn_pcolor_I_surr_out_I2E.pdf")
    plotOutgoing(I_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # Out of curiosity: plot the weights to one neuron (incoming)
    # E-->I
    fname = se.setFName("conn_pcolor_I_surr_in_E2I.pdf")
    plotIncoming(I_surrSp, "I", IneuronIdx, fname,
            ylabel='',
            xlabel='')

    # I-->E
    fname = se.setFName("conn_pcolor_I_surr_in_I2E.pdf")
    plotIncoming(I_surrSp, "E", EneuronIdx, fname,
            ylabel='',
            xlabel='')



