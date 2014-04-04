import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, pcolormesh, subplot2grid, savefig,\
        colorbar, axis, xlabel, ylabel
import matplotlib.ticker as ti
import os
import errno

from interface       import DictDSVisitor
from data_storage.sim_models.ei import extractStateVariable, sumAllVariables, \
        MonitoredSpikes
from plotting.signal import signalPlot
from analysis.spikes import PopulationSpikes
from analysis.grid_cells import extractSpikePositions2D, SNSpatialRate2D, \
        SNAutoCorr, cellGridnessScore
from plotting.bumps  import torusFiringRate
from plotting.grids  import plotSpikes2D
from otherpkg.log    import log_warn, log_info

__all__ = ['DetailedPlotVisitor', 'GridPlotVisitor', 'ISIPlotVisitor']


