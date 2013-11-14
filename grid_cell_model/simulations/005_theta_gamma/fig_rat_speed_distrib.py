#
#   fig_rat_speed_distrib.py
#
#   Rodent speed distribution.
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

from scipy.io           import loadmat
from matplotlib.pyplot  import *
from tables             import *

nBins = 40


data = loadmat('../../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat')
pos_x  = data['pos_x'].ravel()
pos_y  = data['pos_y'].ravel()
rat_dt = data['dt'][0][0]

vel_x = np.diff(pos_x)/rat_dt
vel_y = np.diff(pos_y)/rat_dt
spd = np.sqrt(vel_x**2 + vel_y**2)

h = hist(spd, nBins, normed=True)


h5file = openFile('rat_speed_histogram.h5', mode = "w")

h5file.createArray(h5file.root, 'hist_values', h[0])
h5file.createArray(h5file.root, 'hist_bins', h[1])
h5file.close()


