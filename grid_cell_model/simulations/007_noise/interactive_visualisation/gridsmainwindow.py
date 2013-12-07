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
from PySide.QtCore import SIGNAL
from PySide.QtCore import QObject
from PySide.QtGui import QMainWindow, QFileDialog

from ui_gridsmainwindow import Ui_GridsMainWindow


defaultDir = './output_local/even_spacing/grids/0pA'

class GridsMainWindow(QMainWindow, Ui_GridsMainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent) 
        self.setupUi(self)

        # Connections
        self.loadButton.clicked.connect(self.load_dir)
        self.sweepWidget.dataRenewed.connect(self.gridFieldWidget.updateData)
        self.sweepWidget.dataRenewed.connect(self.aCorrWidget.updateData)
        self.gESpinBox.valueChanged.connect(self.gridFieldWidget.changeRow)
        self.gESpinBox.valueChanged.connect(self.aCorrWidget.changeRow)
        self.gISpinBox.valueChanged.connect(self.gridFieldWidget.changeCol)
        self.gISpinBox.valueChanged.connect(self.aCorrWidget.changeCol)

        # Default states of widgets
        self.loadLineEdit.setText(defaultDir)
        self.gridFieldWidget.changeRow(self.gESpinBox.value())
        self.gridFieldWidget.changeCol(self.gISpinBox.value())
        self.aCorrWidget.changeRow(self.gESpinBox.value())
        self.aCorrWidget.changeCol(self.gISpinBox.value())

        # Internal states
        self.inputDir = None


    def load_dir(self):
        '''
        Opens a dialog box and selects a directory
        '''
        startDir = self.loadLineEdit.text()
        directory = QFileDialog.getExistingDirectory(dir=startDir)
        if (directory):
            self.inputDir = directory
            self.loadLineEdit.setText(self.inputDir)
            shape = (self.ySizeSpinBox.value(),  self.xSizeSpinBox.value())
            self.sweepWidget.setDirectory(self.inputDir, shape)


