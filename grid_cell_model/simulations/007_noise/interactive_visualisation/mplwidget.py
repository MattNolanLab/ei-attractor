import numpy as np
from PySide import QtGui, QtCore

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')
rc('font', size=11)

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from parameters.param_space import JobTrialSpace2D
from data_storage.sim_models import ei

from EI_plotting import sweeps, examples
from EI_plotting import aggregate as aggr

class  MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        self.ax = None
        FigureCanvas.__init__(self, self.fig)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        self.updateGeometry()



class MplWidget(QtGui.QWidget):

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)



class SweepWidget(MplWidget):
    # Signals
    dataRenewed = QtCore.Signal(JobTrialSpace2D)
    positionPicked = QtCore.Signal(int, int)

    # Some default settings
    sweepLeft   = 0.15
    sweepBottom = 0.2
    sweepRight  = 0.9
    sweepTop    = 0.9
    cbar_kw= {
        'label'      : None,
        'location'   : 'right',
        'shrink'     : 0.8,
        'pad'        : 0.05}


    def __init__(self, parent=None):
        MplWidget.__init__(self, parent)
        self.canvas.fig.canvas.mpl_connect('pick_event', self.onpick)

        # Sweeps handling internals
        self.ps = None
        self.varList = None
        self.iterList  = ['g_AMPA_total', 'g_GABA_total']
        self.gEMax, self.gIMax = None, None
        self.pickAnnotation = None
        self.r = None
        self.c = None
        self.ax = None


    def setDirectory(self, rootPath, shape):
        self.dataSpace = JobTrialSpace2D(shape, rootPath)
        gE, gI = aggr.computeYX(self.dataSpace, self.iterList,
                normalize=True)
        self.gEMax = gE[-1, 0]
        self.gIMax = gI[0, -1]
        self.dx = self.gIMax / (self.dataSpace.shape[1] - 1)
        self.dy = self.gEMax / (self.dataSpace.shape[0] - 1)
        data = self.dataSpace[0][0][0].data
        self.noise_sigma = ei.getOption(data, 'noise_sigma')

    def plotPickAnnotation(self, r, c):
        ax = self.ax
        if (ax is not None):
            if (len(ax.lines) != 0):
                del ax.lines[0]
            x = (c + 0.5)*self.dx
            y = (r + 0.5)*self.dy
            self.pickAnnotation, = ax.plot([x], [y], 'o', color='black',
                    markerfacecolor='black', markeredgecolor='white', zorder=2)
            self.canvas.draw()

    def onpick(self, event):
        me = event.mouseevent # should be collection
        xd =  me.xdata
        yd = me.ydata
        if (xd is not None and yd is not None):
            print yd, xd
            r = int(yd / self.gEMax * (self.dataSpace.shape[0] - 1))
            c = int(xd / self.gIMax * (self.dataSpace.shape[1] - 1))
            print r, c
            self.setPickPosition(r, c)
            self.positionPicked.emit(r, c)

    @QtCore.Slot(int, int)
    def setPickPosition(self, r, c):
        self.plotPickAnnotation(r, c)
        if (r != self.r or c != self.c):
            self.r = r
            self.c = c


  

class GridSweepWidget(SweepWidget):

    def __init__(self, parent=None):
        SweepWidget.__init__(self, parent)
        self.varList = ['gridnessScore']


    def setDirectory(self, rootPath, shape):
        super(GridSweepWidget, self).setDirectory(rootPath, shape)

        self.canvas.fig.clear()
        self.ax = self.canvas.fig.add_axes(
                Bbox.from_extents(self.sweepLeft, self.sweepBottom, self.sweepRight,
                                  self.sweepTop))
        self.cbar_kw.update({
            'label'      : 'Gridness score',
            'ticks'      : ti.MultipleLocator(0.5)})

        sweeps.plotGridTrial(self.dataSpace, self.varList, self.iterList,
                self.noise_sigma,
                trialNumList=[],
                sigmaTitle=False,
                ignoreNaNs=True,
                cbar=True, cbar_kw=self.cbar_kw,
                ax=self.ax,
                picker=True)

        c = self.ax.collections
        if (len(c) != 1):
            raise RuntimeError("Something went wrong! len(c) != 1")
        self.dataCollection = c[0]

        self.canvas.draw()
        self.dataRenewed.emit(self.dataSpace)


class BumpSweepWidget(SweepWidget):
    bumpSigmaUpdate = QtCore.Signal(str)
    NTrials = 5
    bumpTStart = 0.5e3

    def __init__(self, parent=None):
        SweepWidget.__init__(self, parent)
        self.types = ['bump_full', 'sigma']
        self.trial = None


    def setDirectory(self, rootPath, shape):
        super(BumpSweepWidget, self).setDirectory(rootPath, shape)

        sigmaBumpText = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'
        self.cbar_kw.update(dict(
                label       = sigmaBumpText,
                ticks       = ti.MultipleLocator(0.1)))

        self.canvas.fig.clear()
        self.ax = self.canvas.fig.add_axes(
                Bbox.from_extents(self.sweepLeft, self.sweepBottom, self.sweepRight,
                                  self.sweepTop))
       
        self.aggrData = aggr.AggregateBumpReciprocal(self.dataSpace, self.iterList,
                self.NTrials, tStart=self.bumpTStart)
        sweeps.plotSweep(self.aggrData,
                self.noise_sigma,
                sigmaTitle=False,
                cbar=True, cbar_kw=self.cbar_kw,
                ax=self.ax,
                picker=True)

        c = self.ax.collections
        if (len(c) != 1):
            raise RuntimeError("Something went wrong! len(c) != 1")
        self.dataCollection = c[0]

        self.canvas.draw()
        self.dataRenewed.emit(self.dataSpace)


    @QtCore.Slot(int)
    def setTrial(self, newTrial):
        self.trial = newTrial
        self.setPickPosition(self.r, self.c)

    @QtCore.Slot(int, int)
    def setPickPosition(self, r, c):
        super(BumpSweepWidget, self).setPickPosition(r, c)
        if (r is not None and c is not None and self.aggrData is not None):
            trialData = self.aggrData.getTrialData()
            print "emitting..."
            self.bumpSigmaUpdate.emit("{0:.3f}".format(trialData[r, c,
                self.trial]))


##############################################################################
# Examples

class ExampleWidget(MplWidget):

    def __init__(self, parent=None):
        MplWidget.__init__(self, parent)
        self.dataSpace = None
        self.r = None
        self.c = None
        self.trial = None
        self.iterList  = ['g_AMPA_total', 'g_GABA_total']

    def clear(self):
        self.canvas.fig.clear()
    

    ########################################################################
    @QtCore.Slot(JobTrialSpace2D)
    def updateData(self, dataSpace):
        self.dataSpace = dataSpace
        self.r = None
        self.c = None


    @QtCore.Slot(int, int)
    def changeRC(self, r, c):
        if (r != self.r or c != self.c):
            self.r = r
            self.c = c
            self.update()

    @QtCore.Slot(int)
    def changeTrial(self, newTrial):
        if (self.trial != newTrial):
            self.trial = newTrial
            self.update()



class GridFieldWidget(ExampleWidget):
    def __init__(self, parent=None):
        ExampleWidget.__init__(self, parent)


    def update(self):
        self.clear()
        rc = (self.r, self.c)
        gsCoords = (0.01, 0.01, 0.99, 0.85)
        if self.dataSpace is not None and rc < tuple(self.dataSpace.shape):
            gs = examples.plotOneGridExample(self.dataSpace, rc, self.iterList,
                    trialNum=self.trial,
                    fig=self.canvas.fig,
                    gsCoords=gsCoords,
                    xlabel=False, ylabel=False,
                    xlabel2=False, ylabel2=False)
        self.canvas.draw()
       

class ACorrelationWidget(ExampleWidget):
    def __init__(self, parent=None):
        ExampleWidget.__init__(self, parent)


    def update(self):
        self.clear()
        rc = (self.r, self.c)
        if self.dataSpace is not None and rc < tuple(self.dataSpace.shape):
            ax = self.canvas.fig.add_subplot(111)
            examples.plotOneGridACorrExample(self.dataSpace, rc, 
                    trialNum=self.trial,
                    ax=ax)
        self.canvas.fig.subplots_adjust(bottom=0, left=0, right=0.99, top=0.85)
        self.canvas.draw()
       

class PopulationFRWidget(ExampleWidget):

    def __init__(self, types, parent=None):
        self.types = types
        ExampleWidget.__init__(self, parent)


    def update(self):
        self.clear()
        rc = (self.r, self.c)
        gsCoords = (0, 0, 1, 0.90)
        if (self.dataSpace is not None and rc < tuple(self.dataSpace.shape)):
            examples.plotOneBumpExample(self.dataSpace, rc, self.iterList,
                    self.types,
                    trialNum=self.trial,
                    exGsCoords=gsCoords,
                    fig=self.canvas.fig)
        self.canvas.draw()


class EFRWidget(PopulationFRWidget):
    def __init__(self, parent=None):
        PopulationFRWidget.__init__(self, ['bump_full', 'rateMap_e'], parent)

class IFRWidget(PopulationFRWidget):
    def __init__(self, parent=None):
        PopulationFRWidget.__init__(self, ['bump_full', 'rateMap_i'], parent)

