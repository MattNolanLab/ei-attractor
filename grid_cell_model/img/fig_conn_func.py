#
#   fig_conn_func.py
#
#   Print connection function
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

from matplotlib.pyplot  import *
from tables             import *


dx = 0.001
x0 = 0
x1 = 1
d = np.arange(x0, x1, dx)

y_dim = np.sqrt(3)/2.
pAMPA_mu = y_dim/2.
pAMPA_sigma = 0.5/6
pGABA_sigma = 0.5/6
pGABA_const = 0.1

shift = 0.1


exc_profile         = np.exp(-(d - pAMPA_mu)**2/2/pAMPA_sigma**2)
exc_profile_shifted = np.exp(-(d - pAMPA_mu - shift)**2/2/pAMPA_sigma**2)
inh_profile         = (1-pGABA_const)*np.exp(-d**2/2./pGABA_sigma**2) + pGABA_const

rcParams['font.size'] = 18

figure(figsize=(6, 4))
plot(d, exc_profile, linewidth=3)
hold('on')
plot(d, inh_profile, linewidth=3)
xlabel('Distance on torus (norm.)')
ylabel('Conn strength (norm.)')
legend(('Excitatory', 'Inhibitory'))
savefig('conn_weights.png')

h5file = openFile('conn_weights_func.h5', mode = "w")

h5file.createArray(h5file.root, 'd', d)
h5file.createArray(h5file.root, 'exc_profile', exc_profile)
h5file.createArray(h5file.root, 'exc_profile_shifted', exc_profile_shifted)
h5file.createArray(h5file.root, 'inh_profile', inh_profile)

h5file.close()

