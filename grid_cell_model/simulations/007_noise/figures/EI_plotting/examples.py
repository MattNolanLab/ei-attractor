#
#   examples.py
#
#   Plotting of E/I examples.
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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from matplotlib.patches  import Rectangle
from matplotlib.gridspec import GridSpec


from . import xlabelText, ylabelText
from . import aggregate as aggr
from plotting.global_defs import globalAxesSettings
from plotting.grids       import plotGridRateMap, plotAutoCorrelation
from plotting.low_level   import xScaleBar
from plotting.bumps       import plotBump
from data_storage.sim_models.ei import extractSummedSignals



def plotOneGridExample(dataSpace, rc, iterList, **kw):
    r, c = rc[0], rc[1]
    spaceRect = (c, r, c, r)
    gsCoords  = kw.pop('gsCoords', (0, 0, 1, 1))
    return drawGridExamples(dataSpace, spaceRect, iterList, gsCoords, **kw)


def plotOneGridACorrExample(dataSpace, rc, trialNum=0, **kw):
    ax = kw.pop('ax', plt.gca())
    kw['rasterized'] = True

    r, c, = rc[0], rc[1]
    trialNumList = None

    corr = dataSpace.aggregateData(['analysis', 'corr'],
            trialNumList=trialNumList, saveData=False, loadData=True,
            output_dtype='list')
    corr_X = dataSpace.aggregateData(['analysis', 'corr_X'], trialNumList=[0],
            saveData=False, loadData=True, output_dtype='list')
    corr_Y = dataSpace.aggregateData(['analysis', 'corr_Y'], trialNumList=[0],
            saveData=False, loadData=True, output_dtype='list')
    arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
            trialNumList=[0], saveData=False, loadData=True,
            output_dtype='array')

    corr      = corr[r][c][trialNum]
    corr_X    = corr_X[r][c][0]
    corr_Y    = corr_Y[r][c][0]
    arenaDiam = arenaDiams[r][c][0]

    plotAutoCorrelation(corr, corr_X, corr_Y, diam=arenaDiam, ax=ax, **kw)


def plotOneBumpExample(sp, rc, iterList, types, **kw):
    #keyword
    wspace    = kw.pop('wspace', 0)
    hspace    = kw.pop('hspace', 0)
    gsCoords  = kw.pop('exGsCoords', (0, 0, 1, 1))

    r, c = rc[0], rc[1] 
    spaceRect = [c, r, c, r]
    return drawBumpExamples(sp, spaceRect, iterList, gsCoords, types,
            xlabel=False, ylabel=False,
            xlabel2=False, ylabel2=False,
            fontsize='x-small',
            rateYPos=1.05, rateXPos=0.98,
            **kw)



##############################################################################

def drawGridExamples(dataSpace, spaceRect, iterList, gsCoords, trialNum=0,
        exIdx=(0, 0), xlabel=True, ylabel=True, xlabelPos=-0.2, xlabel2=True,
        ylabel2=True, ylabelPos=-0.2, xlabel2Pos=-0.6, ylabel2Pos=-0.6,
        fontSize=None, maxRate=True, plotGScore=True, fig=plt.gcf()):
    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = spaceRect[2]
    top    = spaceRect[3]
    exRow, exCol = exIdx

    rateMaps   = aggr.aggregate2D(dataSpace, ['rateMap_e'])
    rateMaps_X = aggr.aggregate2D(dataSpace, ['rateMap_e_X'])
    rateMaps_Y = aggr.aggregate2D(dataSpace, ['rateMap_e_Y'])
    arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
            trialNumList=[0], saveData=False, loadData=True,
            output_dtype='array')
    G          = aggr.aggregate2D(dataSpace, ['gridnessScore'])

    scaleBar = None
    exRows = top - bottom + 1
    exCols = right - left + 1
    gs = GridSpec(exRows, exCols)
    gsLeft   = gsCoords[0]
    gsBottom = gsCoords[1]
    gsRight  = gsCoords[2]
    gsTop    = gsCoords[3]
    gs.update(left=gsLeft, bottom=gsBottom, right=gsRight, top=gsTop)

    we, wi = aggr.computeYX(dataSpace, iterList, r=exRow, c=exCol)
    ax = None
    for r in range(bottom, top+1):
        for c in range(left, right+1):
            print r, c
            rateMap   = rateMaps[r][c][trialNum]
            if (not isinstance(rateMap, np.ndarray)):
                continue

            gsRow = top - r
            gsCol = c - left
            ax = fig.add_subplot(gs[gsRow, gsCol]) 
            X         = rateMaps_X[r][c][0]
            Y         = rateMaps_Y[r][c][0]
            arenaDiam = arenaDiams[r][c][0]
            if (plotGScore):
                gScore = G[r][c][0]
            else:
                gScore = None
            plotGridRateMap(rateMap, X, Y, diam=arenaDiam, scaleBar=scaleBar,
                    scaleText=False, maxRate=maxRate, G=gScore,
                    rasterized=True, ax=ax)

            if (ylabel and gsCol == 0):
                label = "{0:.2f}".format(we[r][c])
                ax.text(ylabelPos, 0.5, label, rotation=90,
                        transform=ax.transAxes, va='center', ha='right',
                        fontsize=fontSize)
            if (xlabel and gsRow == exRows - 1):
                label = "{0:.2f}".format(wi[r][c])
                ax.text(0.5, xlabelPos, label, transform=ax.transAxes,
                        va='top', ha='center', fontsize=fontSize)
            # Second Y label
            if (ylabel2 and r-bottom == 0 and c-left == 0):
                trans = transforms.blended_transform_factory(ax.transAxes,
                        plt.gcf().transFigure)
                weTxt_y = gsBottom + (gsTop - gsBottom)/2.0
                ax.text(ylabel2Pos, weTxt_y, ylabelText, transform=trans,
                        va='center', ha='right', rotation=90,
                        fontsize=fontSize)
            

        # Second X label
        if (xlabel2 and r - bottom == 0):
            trans = transforms.blended_transform_factory(plt.gcf().transFigure,
                    ax.transAxes)
            wiTxt_x = gsLeft + (gsRight - gsLeft)/2.0
            ax.text(wiTxt_x, xlabel2Pos, xlabelText, transform=trans, va='top',
                    ha='center', fontsize=fontSize)

    return gs


def drawBumpExamples(dataSpace, spaceRect, iterList, gsCoords, types, **kw):
    '''
    TODO: code duplication
    '''
    # kw processing
    trialNum   = kw.pop('trialNum', 0)
    exIdx      = kw.pop('exIdx', (0, 0))
    xlabel     = kw.pop('xlabel', True)
    ylabel     = kw.pop('ylabel', True)
    xlabelPos  = kw.pop('xlabelPos', -0.2)
    ylabelPos  = kw.pop('ylabelPos', -0.2)
    xlabel2    = kw.pop('xlabel2', True)
    ylabel2    = kw.pop('ylabel2', True)
    xlabel2Pos = kw.pop('xlabel2os', -0.6)
    ylabel2Pos = kw.pop('ylabel2os', -0.6)
    fontSize   = kw.pop('fontSize', 'medium')
    fig        = kw.pop('fig', plt.gcf())
    # + plotBump() kwargs
    kw['rasterized'] = kw.pop('rasterized', True)
    kw['vmin']       = kw.get('vmin', 0)


    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = spaceRect[2]
    top    = spaceRect[3]
    exRow, exCol = exIdx

    bumps, wi, we = aggr.aggregateType(dataSpace, iterList, types, trialNum+1,
            ignoreNaNs=False, normalizeTicks=False)

    scaleBar = None
    exRows = top - bottom + 1
    exCols = right - left + 1
    gs = GridSpec(exRows, exCols)
    gsLeft   = gsCoords[0]
    gsBottom = gsCoords[1]
    gsRight  = gsCoords[2]
    gsTop    = gsCoords[3]
    gs.update(left=gsLeft, bottom=gsBottom, right=gsRight, top=gsTop)

    ax = None
    for r in range(bottom, top+1):
        for c in range(left, right+1):
            print r, c
            bump = bumps[r][c][trialNum]
            if (not isinstance(bump, np.ndarray)):
                continue

            gsRow = top - r
            gsCol = c - left
            ax = fig.add_subplot(gs[gsRow, gsCol]) 
            plotBump(ax, bump, **kw)

            if (ylabel and gsCol == 0):
                label = "{0:.2f}".format(we[r][c])
                ax.text(ylabelPos, 0.5, label, rotation=90,
                        transform=ax.transAxes, va='center', ha='right',
                        fontsize=fontSize)
            if (xlabel and gsRow == exRows - 1):
                label = "{0:.2f}".format(wi[r][c])
                ax.text(0.5, xlabelPos, label, transform=ax.transAxes,
                        va='top', ha='center', fontsize=fontSize)
            # Second Y label
            if (ylabel2 and r-bottom == 0 and c-left == 0):
                trans = transforms.blended_transform_factory(ax.transAxes,
                        plt.gcf().transFigure)
                weTxt_y = gsBottom + (gsTop - gsBottom)/2.0
                ax.text(ylabel2Pos, weTxt_y, ylabelText, transform=trans,
                        va='center', ha='right', rotation=90,
                        fontsize=fontSize)
            

        # Second X label
        if (xlabel2 and r - bottom == 0):
            trans = transforms.blended_transform_factory(plt.gcf().transFigure,
                    ax.transAxes)
            wiTxt_x = gsLeft + (gsRight - gsLeft)/2.0
            ax.text(wiTxt_x, xlabel2Pos, xlabelText, transform=trans, va='top',
                    ha='center', fontsize=fontSize)

    return gs


def plotSquareGridExample(exLeft, exBottom, sz, fileName, exIdx, sweep_ax,
        sweepDataSpace, iterList, exGsCoords, wspace=0, hspace=0,
        figSize=(2.1,2.1), **kw):
    # Create the example plot
    fig = plt.figure(figsize=figSize)

    exRect = [exLeft, exBottom, exLeft+sz-1, exBottom+sz-1]
    gs = drawGridExamples(sweepDataSpace, exRect, iterList,
            gsCoords=exGsCoords, exIdx=exIdx, fontSize='small', **kw)
    gs.update(wspace=wspace, hspace=hspace)
    plt.savefig(fileName, dpi=300, transparent=False)
    plt.close()

    # Draw the selection into the EI plot
    if (sweep_ax is not None):
        exRow, exCol = exIdx
        Y, X = aggr.computeYX(sweepDataSpace, iterList, r=exRow, c=exCol)
        drawEIRectSelection(sweep_ax, exRect, X, Y)



def drawEIRectSelection(ax, spaceRect, X, Y, color='black'):
    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = spaceRect[2]
    top    = spaceRect[3]

    rectLeft   = X[bottom, left]
    rectBottom = Y[bottom, left]
    rectRight  = X[top, right+1]
    rectTop    = Y[top+1, right]
    #If left==0, we need to shrink the rectangle a little
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    if (left == 0):
        rectLeft += 0.2*dx
    ax.add_patch(Rectangle((rectLeft, rectBottom), rectRight-rectLeft,
        rectTop-rectBottom, facecolor='None', lw=1, edgecolor=color))




def plotGammaExample(ps, r, c, trialNum, tStart, tEnd, **kw):
    ax             = kw.pop('ax', plt.gca())
    noise_sigma    = kw.pop('noise_sigma', None)
    noise_sigma_xy = kw.pop('noise_sigma_xy', (0.95, 0.85))
    xscale_kw      = kw.pop('xscale_kw', None)
    yscale_kw      = kw.pop('yscale_kw', None)

    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    data = ps[r][c][trialNum].data
    monEList = data['stateMon_e']
    t, Isyn_e = extractSummedSignals(monEList, ['I_clamp_GABA_A'], tStart,
            tEnd, monIdx=1)
    monIList = data['stateMon_i']
    t, Isyn_i = extractSummedSignals(monIList, ['I_clamp_AMPA',
        'I_clamp_NMDA'], tStart, tEnd)
    plt.plot(t, Isyn_e, color='red')
    plt.plot(t, Isyn_i, color='blue')
    max = np.max(Isyn_e)
    min = np.min(Isyn_i)
    absmax = np.max([np.abs(max), np.abs(min)])

    ax.set_xlim([tStart, tEnd])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylim([-absmax, absmax])

    if (noise_sigma is not None):
        txt = '$\sigma$ = {0} pA'.format(noise_sigma)
        ns_x, ns_y = noise_sigma_xy[0], noise_sigma_xy[1]
        ax.text(ns_x, ns_y, txt, transform=ax.transAxes, va='bottom', ha='right',
                size='x-small')

    # Scale bars
    if (xscale_kw is not None):
        xscale_kw.update(ax=ax)
        xScaleBar(**xscale_kw)






##############################################################################
# Snapshots of bumps in time

def plotBumpSnapshots(FR, FRt, nSnapshots, **kw):
    fig             = kw.pop('fig')
    timeTitles      = kw.pop('timeTitles', True)
    axesCoords      = kw.pop('axesCoords', (0.12, 0.01, 0.92, 0.7))
    axesDiv         = kw.pop('axesDiv', .01)
    bumpQuality     = kw.pop('bumpQuality', False)
    bumpQualityText = kw.pop('bumpQualityText', '')

    left, bottom, right, top = axesCoords
    width  = right - left
    height = top - bottom

    step = int(FRt.shape[0] / nSnapshots)

    oneWidth = float(width) / nSnapshots
    l = left
    bot = bottom
    indexes = range(0, FRt.shape[0], step)
    max = np.max(FR[:, :, indexes])
    for it in indexes:
        print it
        t = bot + height
        r = l + oneWidth - axesDiv
        print l, bot, r, top

        ax = fig.add_axes(Bbox.from_extents(l, bot, r, top))
        plotBump(ax, FR[:, :, it], vmin=0, vmax=max, rasterized=True, **kw)

        if bumpQuality and it == 0:
            txt = '{0:.2f}'.format(bumpQuality)
            ax.text(-.9, .5, txt, va='center', ha='center',
                    transform=ax.transAxes)

        if timeTitles:
            yTitle = 1.02
            ax.text(.5, yTitle, "{0}".format(FRt[it]*1e-3), size='medium',
                    transform=ax.transAxes, va='bottom', ha='center')
            if it == 0:
                ax.text(.5, yTitle + .3, "t(s)", ha='center', va='bottom',
                        transform=ax.transAxes)
                ax.text(-.9, yTitle, bumpQualityText, ha='center', va='bottom',
                        transform=ax.transAxes)
            #if it / step == nSnapshots - 1:
            #    ax.text(1, yTitle, "s", va='bottom', transform=ax.transAxes)


        l += oneWidth
