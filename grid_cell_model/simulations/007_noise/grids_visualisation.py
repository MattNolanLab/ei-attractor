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
import sys
sys.path.append('/Users/work/work/GridCells/grid_cell_model/simulations/007_noise/figures')

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from PySide.QtGui import QApplication, QFont

from interactive_visualisation.gridsmainwindow import GridsMainWindow


if __name__ == '__main__':
    if sys.platform.startswith('darwin'):
        QFont.insertSubstitution(".Lucida Grande UI", "Lucida Grande")
    app = QApplication(sys.argv)
    window = GridsMainWindow()
    window.show()
    sys.exit(app.exec_())
