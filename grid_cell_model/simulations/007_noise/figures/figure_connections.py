#!/usr/bin/env python
'''
Noise publication figures: connection weight figures.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker     import MultipleLocator, AutoMinorLocator, \
        MaxNLocator
import matplotlib.ticker as ti
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
parser.add_flag('--exampleHists')
parser.add_flag('-w', '--weights')
args = parser.parse_args()

##############################################################################

def plotEToI(sp, gIdx, neuronIdx, trialNum=0, **kw):
    title = kw.pop('title', 'I cell')
    ylim  = kw.pop('ylim', None)

    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_IE']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
            title=title, **kw)
    annG = gE[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_E$ = {0} nS'.format(annG)
    ax.text(0.95, 0.9, ann, ha='right', va='bottom', fontsize='x-small',
            transform=ax.transAxes)
    ax.set_xlim([0, annG])
    ax.set_xticks([0, annG])
    ax.set_ylim(ylim)

def plotIToE(sp, gIdx, neuronIdx, trialNum=0, **kw):
    title = kw.pop('title', 'E cell')

    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]
    ax = pconn.plotConnHistogram(conns,
            title=title, **kw)
    annG = gI[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_I$ = {0} nS'.format(annG)
    ax.text(0.95, 0.9, ann, ha='right', va='bottom', fontsize='x-small',
            transform=ax.transAxes)
    ax.set_xlim([0, annG])
    ax.set_xticks([0, annG])

def plotIToEBrokenAxis(sp, gIdx, neuronIdx, trialNum=0, axBoundaries=[0, 0, 1, 1],
        axesProportions=(0.5, 0.5), bottomLimits=None, topLimits=None,
        **kw):
    title = kw.pop('title', 'E cell')
    fig   = kw.pop('fig', plt.gcf())
    left, bottom, right, top = axBoundaries
    h = top - bottom
    w = right - left
    hBottom = h*axesProportions[0]
    hTop = h*axesProportions[1]
    
    axBottom = fig.add_axes(Bbox.from_extents(left, bottom, right, bottom +
        hBottom))
    axTop = fig.add_axes(Bbox.from_extents(left, top - hTop, right, top),
            sharex=axBottom)

    gE, gI = aggr.computeYX(sp, iterList)
    M      = sp[0][gIdx][trialNum].data['g_EI']
    conns  = M[neuronIdx, :]

    pconn.plotConnHistogram(conns, title=title, ax=axBottom, **kw)
    kw['ylabel'] = ''
    pconn.plotConnHistogram(conns, title=title, ax=axTop, **kw)
    annG = gI[0, gIdx]
    if (annG - int(annG) == 0):
        annG = int(annG)
    ann = '$g_I$ = {0} nS'.format(annG)
    fig.text(left+0.95*w, bottom+0.9*h, ann, ha='right', va='bottom',
            fontsize='x-small')

    axBottom.set_xlim([0, annG])
    axBottom.set_xticks([0, annG])
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


if args.exampleHists or args.all:
    exampleFigSize = (1.6, 1.4)
    exLeft   = 0.4
    exBottom = 0.32
    exRight  = 0.95
    exTop    = 0.85
    exampleRC = ( (5, 15), (15, 5) )
    for exIdx, example in enumerate(exampleRC):
        kw = dict()
        if exIdx == 1:
            kw['xlabel'] = ''

        fig = plt.figure(figsize=exampleFigSize)
        ax = fig.add_axes(Bbox.from_extents(exLeft, exBottom, exRight, exTop))
        plotEToI(sp, example[0], neuronIdx, ylabel='', title='',
                rwidth=0.8,
                linewidth=0,
                **kw)
        #ax.yaxis.set_visible(False)
        ax.yaxis.set_minor_locator(ti.NullLocator())
        #ax.spines['left'].set_visible(False)
        ax.set_xlabel(ax.xaxis.get_label_text(), labelpad=-5)
        #fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=0)
        fname = outputDir + "/figure_connections_examples_E2I{0}.pdf"
        plt.savefig(fname.format(exIdx), dpi=300, transparent=transparent)
        plt.close()


        fig = plt.figure(figsize=exampleFigSize)
        axBoundaries = (exLeft, exBottom, exRight, exTop)
        axBottom, axTop = plotIToEBrokenAxis(sp, example[1], neuronIdx,
                ylabel='', title='',
                axBoundaries=axBoundaries,
                axesProportions=(0.75, 0.2),
                bottomLimits=(0, 60),
                topLimits=(800, 900),
                rwidth=0.8,
                linewidth=0,
                **kw)
        axBottom.set_xlabel(axBottom.xaxis.get_label_text(), labelpad=-5)
        fig.text(exLeft - 0.27, 0.5*(bottom+top), 'Count',
                rotation=90, ha='center', va='center')
        fname = outputDir + "/figure_connections_examples_I2E{0}.pdf"
        plt.savefig(fname.format(exIdx), dpi=300, transparent=transparent)
        plt.close()


if args.weights or args.all:
    # As a control: plot the weights from one neuron (outgoing)
    # E-->I
    fig = plt.figure(figsize=figSize)
    plotOutgoing(sp, gIdx, "E", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_out_E2I.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()

    # I-->E
    fig = plt.figure(figsize=figSize)
    plotOutgoing(sp, gIdx, "I", neuronIdx)
    fig.tight_layout()
    fname = outputDir + "/figure_connections_pcolor_out_I2E.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


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



