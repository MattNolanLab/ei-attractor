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
from matplotlib.ticker  import MaxNLocator, LinearLocator, AutoMinorLocator

from plotting.global_defs import globalAxesSettings, createColorbar
from parameters import DataSpace


###############################################################################

def aggregate2DTrial(sp, varList, trialNumList, fReduce=np.mean):
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
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=True)
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

def plotACTrial(sp, varList, iterList, trialNumList=[0], r=0, c=0, xlabel="",
        ylabel="", colorBar=True, clBarLabel="", vmin=None, vmax=None,
        title="", clbarNTicks=2, xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin, vmax,
            title, clbarNTicks, xticks, yticks)

def plotBumpSigmaTrial(sp, varList, iterList, thr=np.infty, r=0, c=0,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

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
            

def plotGridTrial(sp, varList, iterList, trialNumList=[0], r=0, c=0, xlabel="",
        ylabel="", colorBar=True, clBarLabel="", vmin=None, vmax=None,
        title="", clbarNTicks=2, xticks=True, yticks=True):
    G = aggregate2DTrial(sp, varList, trialNumList)
    G = ma.MaskedArray(G, mask=np.isnan(G))
    Y, X = computeYX(sp, iterList, r=r, c=c)
    return plot2DTrial(X, Y, G, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)


def plot2DTrial(X, Y, C, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="",
        clbarNTicks=2, xticks=True, yticks=True):

    ax = plt.gca()
    globalAxesSettings(ax)
    ax.minorticks_on()
    plt.pcolor(X, Y, C, vmin=vmin, vmax=vmax)
    if (colorBar):
        if (clbarNTicks == None):
            createColorbar(ax, None, clBarLabel, orientation='horizontal',
                    pad=0.2, shrink=0.9)
        else:
            createColorbar(ax, C, clBarLabel, nticks=clbarNTicks,
                    orientation='horizontal', pad=0.2, shrink=0.8)
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

    return C


