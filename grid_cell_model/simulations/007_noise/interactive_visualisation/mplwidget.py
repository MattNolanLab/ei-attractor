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
    # Slots
    dataRenewed = QtCore.Signal(JobTrialSpace2D)

    def __init__(self, parent=None):
        MplWidget.__init__(self, parent)

        # Sweeps handling internals
        self.ps = None
        self.gridVarList = ['gridnessScore']
        self.iterList  = ['g_AMPA_total', 'g_GABA_total']


    def setDirectory(self, rootPath, shape):
        self.dataSpace = JobTrialSpace2D(shape, rootPath)
        data = self.dataSpace[0][0][0].data
        self.noise_sigma = ei.getOption(data, 'noise_sigma')

        sweepLeft   = 0.1
        sweepBottom = 0.2
        sweepRight  = 0.87
        sweepTop    = 0.9
        cbar_kw= {
            'label'      : 'Gridness score',
            'location'   : 'right',
            'shrink'     : 0.8,
            'pad'        : 0.05,
            'ticks'      : ti.MultipleLocator(0.5)}

        self.canvas.fig.clear()
        self.ax = self.canvas.fig.add_axes(
                Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
                                  sweepTop))

        sweeps.plotGridTrial(self.dataSpace, self.gridVarList, self.iterList,
                self.noise_sigma,
                trialNumList=[],
                sigmaTitle=False,
                ignoreNaNs=True,
                cbar=True, cbar_kw=cbar_kw,
                ax=self.ax)
        self.canvas.draw()
        self.dataRenewed.emit(self.dataSpace)


class ExampleWidget(MplWidget):

    def __init__(self, parent=None):
        MplWidget.__init__(self, parent)
        self.dataSpace = None
        self.r = None
        self.c = None
        self.iterList  = ['g_AMPA_total', 'g_GABA_total']

    
    @QtCore.Slot(JobTrialSpace2D)
    def updateData(self, dataSpace):
        self.dataSpace = dataSpace


    @QtCore.Slot(int)
    def changeRow(self, val):
        if (val != self.r):
            self.r = val
            self.update()

    @QtCore.Slot(int)
    def changeCol(self, val):
        if (val != self.c):
            self.c = val
            self.update()

    def clear(self):
        self.canvas.fig.clear()


class GridFieldWidget(ExampleWidget):
    def __init__(self, parent=None):
        ExampleWidget.__init__(self, parent)


    def update(self):
        self.clear()
        rc = (self.r, self.c)
        if self.dataSpace is not None and rc < tuple(self.dataSpace.shape):
            examples.plotOneGridExample(self.dataSpace, rc, self.iterList,
                    fig=self.canvas.fig,
                    xlabel2=False, ylabel2=False)
        self.canvas.draw()
       

class ACorrelationWidget(ExampleWidget):
    def __init__(self, parent=None):
        ExampleWidget.__init__(self, parent)


    def update(self):
        print("acorr update")
        self.clear()
        rc = (self.r, self.c)
        print(rc)
        if self.dataSpace is not None and rc < tuple(self.dataSpace.shape):
            ax = self.canvas.fig.add_subplot(111)
            examples.plotOneGridACorrExample(self.dataSpace, rc, 
                    ax=ax)
        self.canvas.fig.subplots_adjust(bottom=0, left=0, right=1, top=1)
        self.canvas.draw()
       



