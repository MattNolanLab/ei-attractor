#
#   head_direction.py
#
#   Head direction analysis of the sargolini et al., 2006 paper
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

spd_th = 1 # cm/s

dirName = '../../../central_data_store/data/sargolini_et_al_2006/'
#filename = '11084-03020501_t2c1.mat'
filename = '11084-03020501_t2c2.mat'

data = loadmat(dirName + filename)
t = data['t'].flatten()
dt = t[1] - t[0]
x1 = data['x1'].flatten()
x2 = data['x2'].flatten()
y1 = data['y1'].flatten()
y2 = data['y2'].flatten()

# Head direction
hd_x = (x1 - x2)[0:-1]
hd_y = (y1 - y2)[0:-1]

hd = np.arctan2(hd_x, hd_y)

# Movement direction
md_x1 = np.diff(x1)
md_y1 = np.diff(y1)
md1 = np.arctan2(md_x1, md_y1)
spd1 = np.sqrt((md_x1/dt)**2 + (md_y1/dt)**2)
spd1_th_id = spd1 > spd_th

md_x2 = np.diff(x2)
md_y2 = np.diff(y2)
md2 = np.arctan2(md_x2, md_y2)
spd2 = np.sqrt((md_x2/dt)**2 + (md_y2/dt)**2)
spd2_th_id = spd2 > spd_th

figure(figsize=(10, 20))
subplot(2, 1, 1)
plot(hd[spd1_th_id], md1[spd1_th_id], '.')
xlabel('Head direction (rad)')
ylabel('Movement direction LED1')
axis('equal')
subplot(2, 1, 2)
plot(hd[spd2_th_id], md2[spd2_th_id], '.')
xlabel('Head direction (rad)')
ylabel('Movement direction LED2')
axis('equal')

savefig('hd_md_correlation.pdf')

show()


