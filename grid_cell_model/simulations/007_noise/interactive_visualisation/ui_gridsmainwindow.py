# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gridsmainwindow.ui'
#
# Created: Sat Dec  7 00:35:56 2013
#      by: pyside-uic 0.2.15 running on PySide 1.2.1
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_GridsMainWindow(object):
    def setupUi(self, GridsMainWindow):
        GridsMainWindow.setObjectName("GridsMainWindow")
        GridsMainWindow.resize(795, 637)
        self.MainWindow_2 = QtGui.QWidget(GridsMainWindow)
        self.MainWindow_2.setObjectName("MainWindow_2")
        self.verticalLayout = QtGui.QVBoxLayout(self.MainWindow_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_2 = QtGui.QGroupBox(self.MainWindow_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.loadLineEdit = QtGui.QLineEdit(self.groupBox_2)
        self.loadLineEdit.setReadOnly(True)
        self.loadLineEdit.setObjectName("loadLineEdit")
        self.horizontalLayout_2.addWidget(self.loadLineEdit)
        self.loadButton = QtGui.QPushButton(self.groupBox_2)
        self.loadButton.setObjectName("loadButton")
        self.horizontalLayout_2.addWidget(self.loadButton)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.groupBox = QtGui.QGroupBox(self.MainWindow_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.groupBox)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.sweepWidget = SweepWidget(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sweepWidget.sizePolicy().hasHeightForWidth())
        self.sweepWidget.setSizePolicy(sizePolicy)
        self.sweepWidget.setMinimumSize(QtCore.QSize(325, 250))
        self.sweepWidget.setMaximumSize(QtCore.QSize(325, 250))
        self.sweepWidget.setObjectName("sweepWidget")
        self.horizontalLayout_3.addWidget(self.sweepWidget)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setSizeConstraint(QtGui.QLayout.SetFixedSize)
        self.gridLayout.setObjectName("gridLayout")
        spacerItem = QtGui.QSpacerItem(5, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 2, 1, 1)
        self.xSizeSpinBox = QtGui.QSpinBox(self.groupBox)
        self.xSizeSpinBox.setProperty("value", 31)
        self.xSizeSpinBox.setObjectName("xSizeSpinBox")
        self.gridLayout.addWidget(self.xSizeSpinBox, 0, 4, 1, 1)
        self.ySizeSpinBox = QtGui.QSpinBox(self.groupBox)
        self.ySizeSpinBox.setProperty("value", 31)
        self.ySizeSpinBox.setObjectName("ySizeSpinBox")
        self.gridLayout.addWidget(self.ySizeSpinBox, 1, 4, 1, 1)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 0, 6, 1, 1)
        self.label = QtGui.QLabel(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 0, 8, 1, 1)
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 0, 3, 1, 1)
        self.sweepTypeGroupBox = QtGui.QGroupBox(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sweepTypeGroupBox.sizePolicy().hasHeightForWidth())
        self.sweepTypeGroupBox.setSizePolicy(sizePolicy)
        self.sweepTypeGroupBox.setObjectName("sweepTypeGroupBox")
        self.gridLayout_2 = QtGui.QGridLayout(self.sweepTypeGroupBox)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gammaRadioButton = QtGui.QRadioButton(self.sweepTypeGroupBox)
        self.gammaRadioButton.setCheckable(False)
        self.gammaRadioButton.setObjectName("gammaRadioButton")
        self.gridLayout_2.addWidget(self.gammaRadioButton, 0, 0, 1, 1)
        self.gridsRadioButton = QtGui.QRadioButton(self.sweepTypeGroupBox)
        self.gridsRadioButton.setChecked(True)
        self.gridsRadioButton.setObjectName("gridsRadioButton")
        self.gridLayout_2.addWidget(self.gridsRadioButton, 0, 1, 1, 1)
        self.bumpsRadioButton = QtGui.QRadioButton(self.sweepTypeGroupBox)
        self.bumpsRadioButton.setCheckable(False)
        self.bumpsRadioButton.setObjectName("bumpsRadioButton")
        self.gridLayout_2.addWidget(self.bumpsRadioButton, 1, 0, 1, 1)
        self.gridLayout.addWidget(self.sweepTypeGroupBox, 2, 0, 1, 9)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 3, 1, 1)
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 0, 7, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(5, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 0, 5, 1, 1)
        self.gESpinBox = QtGui.QSpinBox(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gESpinBox.sizePolicy().hasHeightForWidth())
        self.gESpinBox.setSizePolicy(sizePolicy)
        self.gESpinBox.setProperty("value", 5)
        self.gESpinBox.setObjectName("gESpinBox")
        self.gridLayout.addWidget(self.gESpinBox, 0, 1, 1, 1)
        self.gISpinBox = QtGui.QSpinBox(self.groupBox)
        self.gISpinBox.setProperty("value", 15)
        self.gISpinBox.setObjectName("gISpinBox")
        self.gridLayout.addWidget(self.gISpinBox, 1, 1, 1, 1)
        self.horizontalLayout_3.addLayout(self.gridLayout)
        self.verticalLayout.addWidget(self.groupBox)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.gridFieldWidget = GridFieldWidget(self.MainWindow_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gridFieldWidget.sizePolicy().hasHeightForWidth())
        self.gridFieldWidget.setSizePolicy(sizePolicy)
        self.gridFieldWidget.setMinimumSize(QtCore.QSize(200, 200))
        self.gridFieldWidget.setMaximumSize(QtCore.QSize(200, 200))
        self.gridFieldWidget.setObjectName("gridFieldWidget")
        self.horizontalLayout.addWidget(self.gridFieldWidget)
        self.aCorrWidget = ACorrelationWidget(self.MainWindow_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.aCorrWidget.sizePolicy().hasHeightForWidth())
        self.aCorrWidget.setSizePolicy(sizePolicy)
        self.aCorrWidget.setMinimumSize(QtCore.QSize(200, 200))
        self.aCorrWidget.setMaximumSize(QtCore.QSize(200, 200))
        self.aCorrWidget.setObjectName("aCorrWidget")
        self.horizontalLayout.addWidget(self.aCorrWidget)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.horizontalLayout)
        GridsMainWindow.setCentralWidget(self.MainWindow_2)
        self.menubar = QtGui.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0, 0, 795, 22))
        self.menubar.setObjectName("menubar")
        GridsMainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(GridsMainWindow)
        self.statusbar.setObjectName("statusbar")
        GridsMainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(GridsMainWindow)
        QtCore.QMetaObject.connectSlotsByName(GridsMainWindow)

    def retranslateUi(self, GridsMainWindow):
        GridsMainWindow.setWindowTitle(QtGui.QApplication.translate("GridsMainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setTitle(QtGui.QApplication.translate("GridsMainWindow", "Input directory", None, QtGui.QApplication.UnicodeUTF8))
        self.loadButton.setText(QtGui.QApplication.translate("GridsMainWindow", "Select", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("GridsMainWindow", "Sweeps", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("GridsMainWindow", "Noise sigma", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("GridsMainWindow", "gE", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("GridsMainWindow", "gI", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("GridsMainWindow", "X size", None, QtGui.QApplication.UnicodeUTF8))
        self.sweepTypeGroupBox.setTitle(QtGui.QApplication.translate("GridsMainWindow", "Sweep type", None, QtGui.QApplication.UnicodeUTF8))
        self.gammaRadioButton.setText(QtGui.QApplication.translate("GridsMainWindow", "Gamma", None, QtGui.QApplication.UnicodeUTF8))
        self.gridsRadioButton.setText(QtGui.QApplication.translate("GridsMainWindow", "Grids", None, QtGui.QApplication.UnicodeUTF8))
        self.bumpsRadioButton.setText(QtGui.QApplication.translate("GridsMainWindow", "Bumps", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("GridsMainWindow", "Y size", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("GridsMainWindow", "??? pA", None, QtGui.QApplication.UnicodeUTF8))

from mplwidget import ACorrelationWidget, GridFieldWidget, SweepWidget
