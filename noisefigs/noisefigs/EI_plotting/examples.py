'''Plotting of various grid field and autocorrelation examples.

.. currentmodule:: noisefigs.EI_plotting.examples

Functions
---------

.. autosummary::

    plotOneGridExample
    plotOneGridACorrExample
    plotOneCorrAngleExample
    plotOneBumpExample
    drawGridExamples
    drawBumpExamples
    plotSquareGridExample
    drawEIRectSelection
    plotGammaExample
    plotBumpSnapshots
'''
from __future__ import absolute_import, print_function, division

from collections import namedtuple

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
from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.plotting.grids       import plotGridRateMap, plotAutoCorrelation
from grid_cell_model.plotting.low_level   import xScaleBar, yScaleBar
from grid_cell_model.plotting.bumps       import plotBump
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals


CorrData = namedtuple('CorrData', ['corr', 'corr_X', 'corr_Y', 'arenaDiams'])


def _getCorrData(dataSpace, populationType):
    '''Get firing field autocorrelation data.

    Parameters
    ----------
    dataSpace : JobTrialSpace2D
        Parameter space to extract data from.
    populationType : str
        Either ``E`` will extract data from the E population, or ``I`` will
        extract data from the I population.

    Returns
    -------
    data : CorrData
        Named tuple that contains the data.
    '''
    trialNumList = None
    if populationType == 'E':
        corr = dataSpace.aggregateData(['analysis', 'corr'],
                                       trialNumList=trialNumList,
                                       saveData=False, loadData=True,
                                       output_dtype='list')
        corr_X = dataSpace.aggregateData(['analysis', 'corr_X'],
                                         trialNumList=[0], saveData=False,
                                         loadData=True, output_dtype='list')
        corr_Y = dataSpace.aggregateData(['analysis', 'corr_Y'],
                                         trialNumList=[0], saveData=False,
                                         loadData=True, output_dtype='list')
        arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
                                             trialNumList=[0], saveData=False,
                                             loadData=True,
                                             output_dtype='array')
    elif populationType == 'I':
        corr = dataSpace.aggregateData(['analysis', 'i_fields', 'corr_i'],
                                       trialNumList=trialNumList,
                                       saveData=False, loadData=True,
                                       output_dtype='list')
        corr_X = dataSpace.aggregateData(['analysis', 'i_fields', 'corr_X'],
                                         trialNumList=[0], saveData=False,
                                         loadData=True, output_dtype='list')
        corr_Y = dataSpace.aggregateData(['analysis', 'i_fields', 'corr_Y'],
                                         trialNumList=[0], saveData=False,
                                         loadData=True, output_dtype='list')
        arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
                                             trialNumList=[0], saveData=False,
                                             loadData=True,
                                             output_dtype='array')
    else:
        raise ValueError("Invalid population type: %s. Use either 'E' or 'I'.",
                         populationType)

    return CorrData(corr=corr, corr_X=corr_X, corr_Y=corr_Y,
                    arenaDiams=arenaDiams)


def plotOneGridExample(dataSpace, rc, iterList, **kw):
    '''Plot an example of one grid firing field.'''
    r, c = rc[0], rc[1]
    spaceRect = (c, r, c, r)
    gsCoords  = kw.pop('gsCoords', (0, 0, 1, 1))
    return drawGridExamples(dataSpace, spaceRect, iterList, gsCoords, **kw)


def plotOneGridACorrExample(dataSpace, rc, trialNum=0, **kw):
    '''Plot an example of one autocorrelation of a grid firing field.'''
    populationType = kw.pop('populationType', 'E')
    ax = kw.pop('ax', plt.gca())
    kw['rasterized'] = True

    r, c, = rc[0], rc[1]

    data = _getCorrData(dataSpace, populationType)

    corr      = data.corr[r][c][trialNum]
    corr_X    = data.corr_X[r][c][0]
    corr_Y    = data.corr_Y[r][c][0]
    arenaDiam = data.arenaDiams[r][c][0]

    plotAutoCorrelation(corr, corr_X, corr_Y, diam=arenaDiam, ax=ax, **kw)


def _getCorrAngleData(dataSpace, populationType):
    '''Get data containing the correlation values and rotation angles.

    Parameters
    ----------
    dataSpace : JobTrialSpace2D
        Parameter space to extract the data from.
    populationType : str
        Either ``E`` for the E cell population, or ``I`` for I cell population.

    Returns
    -------
    correlations, angles : lists of arrays
        Data for correlations and corresponding angles.
    '''
    dataRoot = ['analysis']

    if populationType == 'E':
        pass
    elif populationType == 'I':
        dataRoot += ['i_fields']
    else:
        raise ValueError("Invalid population type: %s. Must be either 'E' or "
                         "'I'", populationType)

    correlations = dataSpace.getReduction(dataRoot + ['gridnessCorr'])
    angles = dataSpace.getReduction(dataRoot + ['gridnessAngles'])

    return correlations, angles


def plotOneCorrAngleExample(dataSpace, rc, trialNum=0, **kw):
    '''A plot of one example of correlation coefficients of rotated
    autocorrelations.
    '''
    populationType = kw.pop('populationType', 'E')
    ax = kw.pop('ax', plt.gca())

    r, c = rc[0], rc[1]
    print("Corr. angle. example; (r, c): ", r, c)

    globalAxesSettings(ax)
    correlations, angles = _getCorrAngleData(dataSpace, populationType)
    ax.plot(angles[r][c][trialNum], correlations[r][c][trialNum])
    ax.set_xlabel('Rotation angle (degrees)')
    ax.set_ylabel('Correlation')
    ax.xaxis.set_major_locator(ti.MultipleLocator(30))
    ax.xaxis.set_major_locator(ti.MultipleLocator(30))


def plotOneBumpExample(sp, rc, iterList, types, **kw):
    '''An example of a bump snapshot for one item in the parameter space.'''
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
RateMapData = namedtuple('RateMapData', ['rateMaps',
                                         'rateMaps_X',
                                         'rateMaps_Y',
                                         'arenaDiams',
                                         'G'])

def _getRateMaps(dataSpace, populationType):
    '''Return the rate map and gridness score data based on the population
    type.

    Parameters
    ----------
    dataSpace : JobTrialSpace2D
        Data space to extract the data from.
    populationType : str
        Either ``E`` for extracting the data from the E population, or ``I``,
        for extracting the data from the E population.

    Returns
    -------
    data : RateMapData
        The named tuple object containing all necessary data.
    '''
    if populationType == 'E':
        rateMaps   = aggr.aggregate2D(dataSpace, ['rateMap_e'])
        rateMaps_X = aggr.aggregate2D(dataSpace, ['rateMap_e_X'])
        rateMaps_Y = aggr.aggregate2D(dataSpace, ['rateMap_e_Y'])
        arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
                                            trialNumList=[0], saveData=False,
                                            loadData=True, output_dtype='array')
        G = aggr.aggregate2D(dataSpace, ['gridnessScore'])
    elif populationType == 'I':
        rateMaps   = aggr.aggregate2D(dataSpace, ['i_fields', 'rateMap_i'])
        rateMaps_X = aggr.aggregate2D(dataSpace, ['i_fields', 'rateMap_i_X'])
        rateMaps_Y = aggr.aggregate2D(dataSpace, ['i_fields', 'rateMap_i_Y'])
        arenaDiams = dataSpace.aggregateData(['options', 'arenaSize'],
                                            trialNumList=[0], saveData=False,
                                            loadData=True, output_dtype='array')
        G = aggr.aggregate2D(dataSpace, ['i_fields', 'gridnessScore'])
    else:
        raise ValueError("Unknown population type: %s, must be either 'E' or "
                         "'I'", populationType)

    return RateMapData(rateMaps=rateMaps, rateMaps_X=rateMaps_X,
                       rateMaps_Y=rateMaps_Y, arenaDiams=arenaDiams, G=G)


def drawGridExamples(dataSpace, spaceRect, iterList, gsCoords, trialNum=0,
        exIdx=(0, 0), xlabel=True, ylabel=True, xlabelPos=-0.2, xlabel2=True,
        ylabel2=True, ylabelPos=-0.2, xlabel2Pos=-0.6, ylabel2Pos=-0.6,
        fontSize=None, maxRate=True, rateStr='Hz', plotGScore=True,
        rasterized=True, fig=plt.gcf(), populationType='E'):
    '''Draw a grid-like plot of grid field examples.'''
    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = min(spaceRect[2], dataSpace.shape[0] - 1)
    top    = min(spaceRect[3], dataSpace.shape[1] - 1)
    exRow, exCol = exIdx

    data = _getRateMaps(dataSpace, populationType)

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
            print(r, c)

            gsRow = top - r
            gsCol = c - left

            ax = fig.add_subplot(gs[gsRow, gsCol])
            X         = data.rateMaps_X[r][c][0]
            Y         = data.rateMaps_Y[r][c][0]

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

            rateMap   = data.rateMaps[r][c][trialNum]
            if not isinstance(rateMap, np.ndarray) or np.isnan(data.G[r][c][0]):
                ax.axis('off')  # remove empty Axes when data is missing
                continue

            arenaDiam = data.arenaDiams[r][c][0]
            if (plotGScore):
                gScore = data.G[r][c][0]
            else:
                gScore = None
            plotGridRateMap(rateMap, X, Y, diam=arenaDiam, scaleBar=scaleBar,
                            scaleText=False, maxRate=maxRate, rateStr=rateStr,
                            G=gScore, rasterized=rasterized, ax=ax, ann_div=.05)



        # Second X label
        if (xlabel2 and r - bottom == 0):
            trans = transforms.blended_transform_factory(plt.gcf().transFigure,
                    ax.transAxes)
            wiTxt_x = gsLeft + (gsRight - gsLeft)/2.0
            ax.text(wiTxt_x, xlabel2Pos, xlabelText, transform=trans, va='top',
                    ha='center', fontsize=fontSize)

    return gs


def drawBumpExamples(dataSpace, spaceRect, iterList, gsCoords, types, **kw):
    '''Draw a grid-like plot of bump attractor (network activity) examples.

    .. todo::
        code duplication
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
            print(r, c)
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
    '''Draw grid field examples into a square grid.'''
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
    '''Draw a rectangular box on top of a parameter sweep.'''
    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = min(spaceRect[2], X.shape[1] - 1)
    top    = min(spaceRect[3], Y.shape[0] - 1)

    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    rectLeft   = X[bottom, left] - .5 * dx
    rectBottom = Y[bottom, left] - .5 * dy
    rectRight  = X[top, right] + .5 * dx
    rectTop    = Y[top, right] + .5 * dy

    ax.add_patch(Rectangle((rectLeft, rectBottom), rectRight-rectLeft,
        rectTop-rectBottom, facecolor='None', lw=1, edgecolor=color))


def plotGammaExample(ps, r, c, trialNum, tStart, tEnd, **kw):
    '''Plot examples of gamma activity.'''
    ax             = kw.pop('ax', plt.gca())
    noise_sigma    = kw.pop('noise_sigma', None)
    noise_sigma_xy = kw.pop('noise_sigma_xy', (0.95, 0.85))
    xscale_kw      = kw.pop('xscale_kw', None)
    yscale_kw      = kw.pop('yscale_kw', None)
    monIdx_e       = kw.pop('monIdx_e', 1)
    monIdx_i       = kw.pop('monIdx_i', 0)

    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    data = ps[r][c][trialNum].data
    monEList = data['stateMon_e']
    t, Isyn_e = extractSummedSignals(monEList, ['I_clamp_GABA_A'], tStart,
                                     tEnd, monIdx=monIdx_e)
    monIList = data['stateMon_i']
    t, Isyn_i = extractSummedSignals(monIList,
                                     ['I_clamp_AMPA', 'I_clamp_NMDA'],
                                     tStart, tEnd,
                                     monIdx=monIdx_i)
    Isyn_e *= 1e-3
    Isyn_i *= 1e-3
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
    if xscale_kw is not None:
        xscale_kw = dict(xscale_kw)
        xscale_kw.update(ax=ax)
        xScaleBar(**xscale_kw)

    if yscale_kw is not None:
        yscale_kw = dict(yscale_kw)
        yscale_kw.update(ax=ax)
        yScaleBar(**yscale_kw)


def plotBumpSnapshots(FR, FRt, tstep, **kw):
    '''Snapshots of bumps in time.'''
    fig             = kw.pop('fig')
    timeTitles      = kw.pop('timeTitles', True)
    axesCoords      = kw.pop('axesCoords', (0.12, 0.01, 0.92, 0.7))
    axesDiv         = kw.pop('axesDiv', .01)
    bumpQuality     = kw.pop('bumpQuality', False)
    bumpQualityText = kw.pop('bumpQualityText', '')
    bumpQualityX    = kw.pop('bumpQualityX', -.9)
    maxRateColor    = kw.pop('maxRateColor', 'w')

    left, bottom, right, top = axesCoords
    width  = right - left
    height = top - bottom

    indexes = range(0, FRt.shape[0], tstep)
    oneWidth = float(width) / len(indexes)
    l = left
    bot = bottom
    max = np.max(FR[:, :, indexes])
    lastIndex = len(indexes) - 1
    idx = 0
    for it in indexes:
        print(it)
        t = bot + height
        r = l + oneWidth - axesDiv
        print(l, bot, r, top)

        ax = fig.add_axes(Bbox.from_extents(l, bot, r, top))
        plotBump(ax, FR[:, :, it], vmin=0, vmax=max, rasterized=True, **kw)
        if idx == lastIndex:
            rateText = "%.0f Hz" % max
            ax.text(1.05, .95, rateText, ha='left', va='top',
                    color=maxRateColor, transform=ax.transAxes, size='small',
                    weight='bold', clip_on=False)

        if bumpQuality and it == 0:
            txt = '{0:.2f}'.format(bumpQuality)
            ax.text(bumpQualityX, .5, txt, va='center', ha='center',
                    transform=ax.transAxes)

        if timeTitles:
            yTitle = 1.02
            ax.text(.5, yTitle, "{0}".format(FRt[it]*1e-3), size='medium',
                    transform=ax.transAxes, va='bottom', ha='center')
            if it == 0:
                ax.text(.5, yTitle + .3, "t(s)", ha='center', va='bottom',
                        transform=ax.transAxes)
                ax.text(bumpQualityX, yTitle, bumpQualityText, ha='center',
                        va='bottom', transform=ax.transAxes)

        l += oneWidth
        idx += 1

    return max  # Hack, but hopefully ok for now
