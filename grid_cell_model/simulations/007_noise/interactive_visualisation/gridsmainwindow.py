#!/usr/bin/env python
#
#   grids_visualisation.py
#
#   A simple grid data visualisation program.
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
from matplotlib import rcParams
rcParams['backend.qt4'] = 'PySide'

import sys
from PySide import QtCore
from PySide.QtCore import SIGNAL
from PySide.QtCore import QObject
from PySide.QtGui import QMainWindow, QFileDialog

from ui_gridsmainwindow import Ui_GridsMainWindow


class GridsMainWindow(QMainWindow, Ui_GridsMainWindow):
    defaultBaseDir = './output_local/even_spacing'
    gridsSubDir = 'grids'
    bumpsSubDir = 'bumps'


    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent) 
        self.setupUi(self)

        # Connections
        self.gridSweepWidget.dataRenewed.connect(self.gridFieldWidget.updateData)
        self.gridSweepWidget.dataRenewed.connect(self.aCorrWidget.updateData)
        self.gridSweepWidget.positionPicked.connect(self.updateRC)

        self.bumpSweepWidget.dataRenewed.connect(self.eFRExampleWidget.updateData)
        self.bumpSweepWidget.dataRenewed.connect(self.iFRExampleWidget.updateData)
        self.bumpSweepWidget.positionPicked.connect(self.updateRC)
        self.bumpSweepWidget.bumpSigmaUpdate.connect(self.bumpSigmaLabel.setText)

        self.tabWidget.currentChanged.connect(self.updateExamples)
        self.noise_sigmaSpinBox.valueChanged.connect(self.changeNoiseSigma)
        self.gESpinBox.valueChanged.connect(self.gESpinBoxChange)
        self.gISpinBox.valueChanged.connect(self.gISpinBoxChange)

        self.selectButton.clicked.connect(self.select_dir)
        self.loadButton.clicked.connect(self.load_dir)

        self.trialNumSpinBox.valueChanged.connect(self.bumpSweepWidget.setTrial)
        self.trialNumSpinBox.valueChanged.connect(self.updateExamples)

        # Default states of widgets
        self.loadLineEdit.setText(self.defaultBaseDir)
        self.useNoiseSigmaCheckBox.setCheckState(QtCore.Qt.Checked)
        self.r = self.gESpinBox.value()
        self.c = self.gISpinBox.value()
        self.bumpSigmaLabel.setText("???")

        # Initial signal distribution
        self.trialNumSpinBox.valueChanged.emit(self.trialNumSpinBox.value())



    def constructFinalDir(self, dirType):
        baseDir = self.loadLineEdit.text()
        subDir = self.decideSubDir()
        return "{0}/{1}/{2}pA".format(baseDir, dirType, subDir)


    def constructFinalDirs(self):
        class FinalDirs(object):
            pass
        retVal = FinalDirs()
        for dirType in ['gamma_bump', 'grids', 'velocity']:
            finalDir = self.constructFinalDir(dirType)
            print finalDir
            setattr(retVal, dirType, finalDir)
        return retVal


    #########################################################################
    # Slots
    def select_dir(self):
        '''
        Opens a dialog box and selects a directory
        '''
        startDir = self.loadLineEdit.text()
        directory = QFileDialog.getExistingDirectory(dir=startDir)
        if (directory):
            self.loadLineEdit.setText(directory)


    def gESpinBoxChange(self, newVal):
        if (self.r != newVal):
            self.r = self.gESpinBox.value()
            self.updateRC(self.r, self.c)

    def gISpinBoxChange(self, newVal):
        if (self.c != newVal):
            self.c = self.gISpinBox.value()
            self.updateRC(self.r, self.c)
        

    @QtCore.Slot(int, int)
    def updateRC(self, r, c):
        self.r, self.c = r, c
        self.gridSweepWidget.setPickPosition(r, c)
        self.bumpSweepWidget.setPickPosition(r, c)
        self.updateExamples()
        self.gESpinBox.setValue(r)
        self.gISpinBox.setValue(c)


    def updateExamples(self):
        r = self.r
        c = self.c
        trial = self.trialNumSpinBox.value()
        currentTab = self.tabWidget.currentIndex()
        if (currentTab == 0):
            self.gridFieldWidget.changeTrial(trial) # Must be first
            self.aCorrWidget.changeTrial(trial)
            self.gridFieldWidget.changeRC(r, c)
            self.aCorrWidget.changeRC(r, c)
        elif (currentTab == 1):
            self.eFRExampleWidget.changeTrial(trial)
            self.iFRExampleWidget.changeTrial(trial)
            self.eFRExampleWidget.changeRC(r, c)
            self.iFRExampleWidget.changeRC(r, c)
        else:
            raise ValueError('No such tab with index: {0}'.format(currentTab))


    def load_dir(self):
        allDirs = self.constructFinalDirs() 
        shape = (self.ySizeSpinBox.value(),  self.xSizeSpinBox.value())
        self.gridSweepWidget.setDirectory(allDirs.grids, shape)
        self.bumpSweepWidget.setDirectory(allDirs.gamma_bump, shape)
        self.updateRC(self.r, self.c)

    def decideSubDir(self):
        if self.useNoiseSigmaCheckBox.checkState() == QtCore.Qt.Checked:
            return self.noise_sigmaSpinBox.value()
        else:
            return self.subDirLineEdit.text()


    #########################################################################
    @QtCore.Slot(int)
    def changeNoiseSigma(self, ns):
        if (self.autoLoadCheckBox.checkState() == QtCore.Qt.Checked and
                self.useNoiseSigmaCheckBox.checkState() == QtCore.Qt.Checked):
            self.load_dir()
            self.updateRC(self.r, self.c)

