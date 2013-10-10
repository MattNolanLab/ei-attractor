#
#   plotting.py
#
#   Shared plotting functions.
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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.ticker   import MaxNLocator, LinearLocator, AutoMinorLocator, \
        MultipleLocator, NullLocator
from matplotlib.patches  import Rectangle
from matplotlib.gridspec import GridSpec


from plotting.global_defs import globalAxesSettings
from plotting.grids       import plotGridRateMap
from parameters import DataSpace
from figures_shared       import plotBump, createColorbar

###############################################################################

xlabelText = '$g_I$ (nS)'
ylabelText = '$g_E$ (nS)'

###############################################################################

def aggregate2DTrial(sp, varList, trialNumList, fReduce=np.mean,
        ignoreNaNs=False):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space.

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    trialNumList : list of ints
        A list specifying exactly which trials are to be processed.
    fReduce : a function f(data, axis)
        A reduction function.
    ignoreNaNs : bool, optional
        If True, mask the NaN values.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=True)
    if (ignoreNaNs):
        nans = np.isnan(retVar)
        retVar = ma.MaskedArray(retVar, mask=nans)
    return fReduce(retVar, 2)


def aggregate2D(sp, varList, funReduce=None):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space, however the data is retrieved
    from the top-level of the data hierarchy, i.e. sp['analysis']. funReduce is
    applied on the data of each item in the ParamSpace (not necessarily
    trials).

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    funReduce : a function f(data, axis)
        A reduction function.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    return sp.aggregateData(varList, funReduce=funReduce,
            trialNumList='all-at-once', saveData=True)



def computeYX(sp, iterList, r=0, c=0, trialNum=0):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ne')
    Ni = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ni')
    return E/Ne, I/Ni



def plotACTrial(sp, varList, iterList, noise_sigma, trialNumList=[0], **kw):
    #kw arguments
    r          = kw.pop('r', 0)
    c          = kw.pop('c', 0)
    cbar       = kw.pop('cbar', True)
    sigmaTitle = kw.pop('sigmaTitle', True)

    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    print("    max(C): {0}".format(np.max(C)))
    print("    min(C): {0}".format(np.min(C)))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    return C, ax, cax


def plotBumpSigmaTrial(sp, varList, iterList, noise_sigma, trialNumList=[0],
        thr=np.infty, **kw):
    #kw arguments
    r          = kw.pop('r', 0)
    c          = kw.pop('c', 0)
    cbar       = kw.pop('cbar', True)
    kw['cmap'] = kw.get('cmap', 'jet_r')
    sigmaTitle = kw.pop('sigmaTitle', True)

    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    C, ax, cax =  plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    return C, ax, cax



def plotBumpErrTrial(sp, varList, iterList, thr=np.infty, r=0, c=0, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    C = np.sqrt(aggregate2DTrial(sp, varList, trialNumList))
    if mask is None:
        mask = False
    C = ma.MaskedArray(C, mask=np.logical_or(np.logical_or(np.isnan(C), C >
        thr), mask))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotFRTrial(sp, varList, iterList, thr=np.infty, r=0, c=0, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    FR = aggregate2DTrial(sp, varList, trialNumList)
    if mask is None:
        mask = False
    FR = ma.MaskedArray(FR, mask=np.logical_or(FR > thr, mask))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    return plot2DTrial(X, Y, FR, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)
            

def plotGridTrial(sp, varList, iterList, noise_sigma, trialNumList=[0], **kwargs):
    #kw arguments
    r          = kwargs.pop('r', 0)
    c          = kwargs.pop('c', 0)
    nansAs0    = kwargs.pop('nansAs0', False)
    ignoreNaNs = kwargs.pop('ignoreNaNs', False)
    sigmaTitle = kwargs.pop('sigmaTitle', True)
    cbar       = kwargs.pop('cbar', True)

    G = aggregate2DTrial(sp, varList, trialNumList, ignoreNaNs=ignoreNaNs)
    nans = np.isnan(G)
    if (nansAs0):
        G[nans] = 0
    else:
        G = ma.MaskedArray(G, mask=nans)
    Y, X = computeYX(sp, iterList, r=r, c=c)
    G, ax, cax = plot2DTrial(X, Y, G, colorBar=cbar, **kwargs)

    print("    max(G): {0}".format(np.max(G)))
    print("    min(G): {0}".format(np.min(G)))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))
    return G, ax, cax


def computeVelYX(sp, iterList, r=0, c=0, trialNum=0):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0], 'net_Ne')
    Ni = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0], 'net_Ni')
    return E/Ne, I/Ni



def plotVelTrial(sp, varList, iterList, noise_sigma, **kwargs):
    # process kwargs
    r          = kwargs.pop('r', 0)
    c          = kwargs.pop('c', 0)
    sigmaTitle = kwargs.pop('sigmaTitle', True)
    cbar       = kwargs.pop('cbar', True)

    C    = np.abs(aggregate2D(sp, varList, funReduce = np.sum))
    C    = ma.MaskedArray(C, mask = np.isnan(C))
    Y, X = computeVelYX(sp, iterList, r=r, c=c)
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kwargs)

    print("plotVelTrial: max(C): {0}".format(np.max(C.ravel())))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    return C, ax, cax


def plotVelStdSweep(sp, iterList, noise_sigma, **kwargs):
    # process kwargs
    r          = kwargs.pop('r', 0)
    c          = kwargs.pop('c', 0)
    sigmaTitle = kwargs.pop('sigmaTitle', True)
    cbar       = kwargs.pop('cbar', True)


    varList = ['analysis', 'bumpVelAll']
    all = sp.aggregateData(varList, funReduce=None,
        trialNumList='all-at-once', output_dtype='list')
    C = np.ndarray((len(all), len(all[0])))
    for row in xrange(len(all)):
        for col in xrange(len(all[0])):
            slopes = all[row][col]
            if (isinstance(slopes, np.ndarray)):
                std = np.std(slopes, axis=0)
                C[row, col] = np.mean(std)
            else:
                C[row, col] = np.nan

    Y, X = computeVelYX(sp, iterList, r=r, c=c)
    C = ma.MaskedArray(C, mask = np.isnan(C))
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kwargs)

    print("plotVelStdSweep: max(C): {0}".format(np.max(C.ravel())))
    print("plotVelStdSweep: min(C): {0}".format(np.min(C.ravel())))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))
    
    return C, ax, cax



def plot2DTrial(X, Y, C, ax=plt.gca(), xlabel=xlabelText, ylabel=ylabelText,
        colorBar=False, clBarLabel="", vmin=None, vmax=None, title="",
        clbarNTicks=2, xticks=True, yticks=True, cmap=None, cbar_kwargs={},
        **kw):
    # kw arguments (cbar)
    cbar_kwargs['label']       = cbar_kwargs.get('label', '')
    cbar_kwargs['shrink']      = cbar_kwargs.get('shrink', 0.8)
    cbar_kwargs['orientation'] = cbar_kwargs.get('orientation', 'vertical')
    cbar_kwargs['pad']         = cbar_kwargs.get('pad', 0.05)
    cbar_kwargs['ticks']       = cbar_kwargs.get('ticks', MultipleLocator(5))

    globalAxesSettings(ax)
    ax.minorticks_on()
    plt.pcolor(X, Y, C, vmin=vmin, vmax=vmax, cmap=cmap, **kw)
    cax = createColorbar(ax, **cbar_kwargs)
    if (colorBar == False):
        cax.set_visible(False)
    if (xlabel != ""):
        plt.xlabel(xlabel, va='top')
        ax.xaxis.set_label_coords(0.5, -0.125)
    if (ylabel != ""):
        plt.ylabel(ylabel, ha='right')
        ax.yaxis.set_label_coords(-0.125, 0.5)
    ax.xaxis.set_ticks([0, 6])
    ax.yaxis.set_ticks([0, 6])
    ax.xaxis.set_minor_locator(AutoMinorLocator(6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(6))
    plt.axis('scaled')
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return C, ax, cax



def drawGridExamples(dataSpace, spaceRect, iterList, gsCoords, trialNum=0,
        exIdx=(0, 0), xlabel=True, ylabel=True, xlabelPos=-0.2, xlabel2=True,
        ylabel2=True, ylabelPos=-0.2, xlabel2Pos=-0.6, ylabel2Pos=-0.6,
        fontSize=None, maxRate=True, plotGScore=True):
    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = spaceRect[2]
    top    = spaceRect[3]
    exRow, exCol = exIdx

    rateMaps   = aggregate2D(dataSpace, ['analysis', 'rateMap_e'])
    rateMaps_X = aggregate2D(dataSpace, ['analysis', 'rateMap_e_X'])
    rateMaps_Y = aggregate2D(dataSpace, ['analysis', 'rateMap_e_Y'])
    arenaDiams = aggregate2D(dataSpace, ['options', 'arenaSize'])
    G          = aggregate2D(dataSpace, ['analysis', 'gridnessScore'])

    scaleBar = None
    exRows = top - bottom + 1
    exCols = right - left + 1
    gs = GridSpec(exRows, exCols)
    gsLeft   = gsCoords[0]
    gsBottom = gsCoords[1]
    gsRight  = gsCoords[2]
    gsTop    = gsCoords[3]
    gs.update(left=gsLeft, bottom=gsBottom, right=gsRight, top=gsTop)

    we, wi = computeYX(dataSpace, iterList, r=exRow, c=exCol)
    ax = None
    for r in range(bottom, top+1):
        for c in range(left, right+1):
            print r, c
            rateMap   = rateMaps[r][c][trialNum]
            if (not isinstance(rateMap, np.ndarray)):
                continue

            gsRow = top - r
            gsCol = c - left
            ax = plt.subplot(gs[gsRow, gsCol]) 
            X         = rateMaps_X[r][c][0]
            Y         = rateMaps_Y[r][c][0]
            arenaDiam = arenaDiams[r][c][0]
            if (plotGScore):
                gScore = G[r][c][0]
            else:
                gScore = None
            plotGridRateMap(rateMap, X, Y, diam=arenaDiam, scaleBar=scaleBar,
                    scaleText=False, maxRate=maxRate, G=gScore)

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


def drawBumpExamples(dataSpace, spaceRect, iterList, gsCoords, **kw):
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
    # + plotBump() kwargs


    left   = spaceRect[0]
    bottom = spaceRect[1]
    right  = spaceRect[2]
    top    = spaceRect[3]
    exRow, exCol = exIdx

    bumps = aggregate2D(dataSpace, ['analysis', 'bump_e', 'bump_e_rateMap'])

    scaleBar = None
    exRows = top - bottom + 1
    exCols = right - left + 1
    gs = GridSpec(exRows, exCols)
    gsLeft   = gsCoords[0]
    gsBottom = gsCoords[1]
    gsRight  = gsCoords[2]
    gsTop    = gsCoords[3]
    gs.update(left=gsLeft, bottom=gsBottom, right=gsRight, top=gsTop)

    we, wi = computeYX(dataSpace, iterList, r=exRow, c=exCol)
    ax = None
    for r in range(bottom, top+1):
        for c in range(left, right+1):
            print r, c
            bump = bumps[r][c][trialNum]
            if (not isinstance(bump, np.ndarray)):
                continue

            gsRow = top - r
            gsCol = c - left
            ax = plt.subplot(gs[gsRow, gsCol]) 
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
        sweepDataSpace, iterList, exGsCoords, wspace=0, hspace=0, figSize=(2.1, 2.1),
        fontSize=None, xlabel=True, ylabel=True, xlabel2=True, ylabel2=True,
        maxRate=True, plotGScore=True):
    # Create the example plot
    fig = plt.figure(figsize=figSize)

    exRect = [exLeft, exBottom, exLeft+sz-1, exBottom+sz-1]
    gs = drawGridExamples(sweepDataSpace, exRect, iterList,
            gsCoords=exGsCoords, exIdx=exIdx, fontSize='small', xlabel=xlabel,
            ylabel=ylabel, xlabel2=xlabel2, ylabel2=ylabel2, maxRate=maxRate,
            plotGScore=plotGScore)
    gs.update(wspace=wspace, hspace=hspace)
    plt.savefig(fileName, dpi=300, transparent=False)
    plt.close()

    # Draw the selection into the EI plot
    if (sweep_ax is not None):
        exRow, exCol = exIdx
        Y, X = computeYX(sweepDataSpace, iterList, r=exRow, c=exCol)
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


##############################################################################
# Slices through parameter spaces
def decideLabels(type):
    if (type == 'horizontal'):
        xlabel = ylabelText
        titleText = "$g_I$"
    elif (type == 'vertical'):
        xlabel = xlabelText
        titleText = "$g_E$"
    else:
        raise ValueError("type must be 'horizontal' or 'vertical'")

    titles = dict(xlabel=xlabel, titleText=titleText)
    return titles


def plotOneSlice(ax, x, y, **kw):
    xlabel           = kw.pop('xlabel', '')
    ylabel           = kw.pop('ylabel', '')
    ylabelPos        = kw.pop('ylabelPos', -0.2)
    fmt              = kw.pop('fmt', 'o-')
    xticks           = kw.pop('xticks', True)
    yticks           = kw.pop('yticks', True)
    errors           = kw.pop('errors', True)
    kw['markersize'] = kw.get('markersize', 4)

    globalAxesSettings(ax)
    ndim = len(y.shape)
    if (ndim == 2):
        mean = np.mean(y, axis=1) # axis 1: trials
        std  = np.std(y, axis=1)
        if (errors):
            ax.errorbar(x, mean, std, fmt=fmt, **kw)
        else:
            ax.plot(x, mean, fmt, **kw)
    elif (ndim == 1):
        ax.plot(x, y, fmt, **kw)
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, rotation=90, transform=ax.transAxes,
            va='center', ha='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    w = x[-1] - x[0]
    margin = 0.025
    ax.set_xlim([-margin*w, x[-1]+margin*w])

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])


def extractSliceX(X, Y, rowSlice, colSlice, type):
    if (type == 'horizontal'):
        return X[0, colSlice]
    else:
        return Y[rowSlice, 0]


def collapseTrials(var, type, ignoreNaNs):
    if (len(var.shape) != 3):
        return var

    trials = []
    for trialIdx in xrange(var.shape[2]):
        trials.append(var[:, :, trialIdx])

    res    = None
    if (type == 'horizontal'):
        res = np.vstack(trials).T
    elif (type == 'vertical'):
        res = np.hstack(trials)
    else:
        raise ValueError("type must be 'horizontal' or 'vertical'")

    if (ignoreNaNs):
        nans = np.isnan(res)
        res = ma.MaskedArray(res, mask=nans)

    return res


def plotGridnessSlice(paramSpaces, rowSlice, colSlice, type, NTrials=1, **kw):
    # kwargs
    title        = kw.pop('title', True)
    rcG          = kw.pop('rowsCols', [(1, 22), (1, 22), (1, 22)]) # (row, col)
    ax           = kw.pop('ax', plt.gca())
    iterList     = kw.pop('iterList', ['g_AMPA_total', 'g_GABA_total'])
    kw['ylabel'] = kw.get('ylabel', 'Gridness score')
    labels = decideLabels(type)
    kw['xlabel'] = kw.get('xlabel', labels['xlabel'])
    ignoreNaNs   = kw.get('ignoreNaNs', True)

    GVars = ['gridnessScore']
    trialNumList = range(NTrials)
    sp = paramSpaces

    # Gridness score
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.grids[idx]
        G = space.aggregateData(GVars, trialNumList, output_dtype='array',
                loadData=True, saveData=False)
        if (ignoreNaNs):
            nans = np.isnan(G)
            G = ma.MaskedArray(G, mask=nans)
        Y, X = computeYX(space, iterList, r=rcG[idx][0], c=rcG[idx][1])
        x = extractSliceX(X, Y, rowSlice, colSlice, type)
        y = collapseTrials(G[rowSlice, colSlice, :], type, ignoreNaNs)
        #import pdb; pdb.set_trace()
        plotOneSlice(ax, x, y, **kw)
    ax.yaxis.set_major_locator(MaxNLocator(4))

    # Title
    if (title):
        if (type == 'horizontal'):
            sliceG = Y[rowSlice, 0]
        else:
            sliceG = X[0, colSlice]
        if (isinstance(sliceG, np.ndarray)):
            txt = '{0} = {1} - {2} nS'.format(labels['titleText'], sliceG[0],
                    sliceG[-1])
        else:
            txt = '{0} = {1} nS'.format(labels['titleText'], sliceG)
        ax.set_title(txt, x=0.99, y=0.92, va='bottom', ha='right',
                fontsize='small')
        
    return ax
        
##############################################################################
# Raster plots
def plotEIRaster(ESpikes, ISpikes, tLimits, **kw):
    # kw arguments 
    ax               = kw.pop('ax', plt.gca())
    ylabel           = kw.pop('ylabel', 'Neuron #')
    yticks           = kw.pop('yticks', True)
    yticks_style     = kw.pop('yticks_style', 'separate')
    ylabelPos        = kw.pop('ylabelPos', -0.2)
    EColor           = kw.pop('ecolor', 'red')
    IColor           = kw.pop('icolor', 'blue')
    kw['markersize'] = kw.get('markersize', 1.0)

    ESpikes = ESpikes.windowed(tLimits)
    ISpikes = ISpikes.windowed(tLimits)

    ESenders, ETimes = ESpikes.rasterData()
    ISenders, ITimes = ISpikes.rasterData()
    ISenders += ESpikes.N

    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_minor_locator(NullLocator())

    ax.plot(ETimes, ESenders+1, '.', color='red',  mec='none', **kw)
    ax.plot(ITimes, ISenders+1, '.', color='blue', mec='none', **kw)

    ax.set_xlim(tLimits)
    ax.set_ylim([1, ESpikes.N+ISpikes.N])
    if (yticks_style == 'separate'):
        ax.set_yticks([1, ESpikes.N, ESpikes.N+ISpikes.N])
    ax.invert_yaxis()
    ax.text(ylabelPos, 0.5, ylabel, va='center', ha='right',
            transform=ax.transAxes, rotation=90)
    if (not yticks):
        ax.yaxis.set_ticklabels([])
    
    return ax
