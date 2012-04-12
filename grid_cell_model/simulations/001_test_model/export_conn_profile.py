#
#   export_conn_profile.py
#
#   Script to export connection profile into HDF5 file. Requires pytables.
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
import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from matplotlib.backends.backend_pdf import PdfPages

from scipy import linspace
from scipy.io import loadmat
from scipy.io import savemat
from optparse import OptionParser
from datetime import datetime

import time
import math
import sys
import numpy as np
import logging as lg

from EI_network import *
from EI_network_sim_mod import *
from custombrian import *

from tools import *
from plotting import *

lg.basicConfig(level=lg.WARNING)


parser = getOptParser()

parser.add_option("--Ivel", type="float", help="Velocity input (pA)")
parser.add_option("--pAMPA_sigma", type="float", help="AMPA profile spread (normalised)")
parser.add_option("--Iext_e_min", type=float,
        help="Minimal external current onto E cells (theta stim.) (A)")
parser.add_option("--Iext_i_min", type=float,
        help="Minimal external current onto I cells (theta stim.) (I)")
parser.add_option("--g_extraGABA_total", type=float,
        help="Uniform inhibition (E-->I only) total conductance (S)")
parser.add_option("--extraGABA_density", type=float,
        help="Uniform inhibition (E-->I only) connection density")
parser.add_option("--prefDirC", type=float,
        help="Preferred directtion multiplier")

(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()


output_fname = 'conn_export_regular_torus.h5'

options.ndim = 2 #'twisted_torus'
ei_net = EI_Network(options, simulationClock)

# Mexican hat properties and AMPA/GABA connections
y_dim = np.sqrt(3)/2.
pAMPA_mu = y_dim/2.
pAMPA_sigma = 0.5/6
pGABA_sigma = 0.5/6
ei_net.connMexicanHat(pAMPA_mu, pAMPA_sigma, pGABA_sigma)
ei_net.randomInhibition(options.g_extraGABA_total, options.extraGABA_density)

print('pAMPA_sigma = ' + str(options.pAMPA_sigma) + '/6')


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################


#ei_net.net.run(0.1*second, report='stdout',
#        report_period=options.update_interval*second)

##outData = {}
###outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
###if spikeMon_e is not None:
###    outData['spikeCell_e'] = spikeMon_e.aspikes
###if spikeMon_i is not None:
###    outData['spikeCell_i'] = spikeMon_i.aspikes
###outData['options'] = options._einet_optdict

##outData['bumpPos'] = pos
##outData['bumpPos_times'] = bumpPos_times

###outData['Fe'] = Fe
###outData['Fe_t'] = Fe_t
###outData['Fi'] = Fi
###outData['Fi_t'] = Fi_t

##outData['theta_stateMon_Iclamp_e_times'] = theta_stateMon_Iclamp_e.times
##outData['theta_stateMon_Iclamp_e_values'] = theta_stateMon_Iclamp_e.values
##outData['theta_stateMon_Iclamp_i_times'] = theta_stateMon_Iclamp_i.times
##outData['theta_stateMon_Iclamp_i_values'] = theta_stateMon_Iclamp_i.values
##
##savemat(output_fname + '_output.mat', outData, do_compression=True)

saveIgor = True
if saveIgor:
    from tables import *
    h5file = openFile(output_fname, mode = "w", title =
            "Attractor export figures")

    if options.ndim == 2:
        ei_net.Ne_x = options.Ne
        ei_net.Ne_y = options.Ne
        ei_net.Ni_x = options.Ni
        ei_net.Ni_y = options.Ni


    # Export connection profiles in the middle of the sheet
    n_e_it = ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1
    n_i_it = ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1
    AMPA_export = np.reshape(ei_net.AMPA_conn.W.todense()[n_e_it, :],
            (ei_net.Ni_y, ei_net.Ni_x))
    GABA_export = np.reshape(ei_net.GABA_conn1.W.todense()[n_i_it, :] +
        ei_net.extraGABA_conn1.W.todense()[n_i_it, :], (ei_net.Ne_y,
            ei_net.Ne_x))

    h5file.createArray(h5file.root, 'AMPA_profile_center', AMPA_export)
    h5file.createArray(h5file.root, 'GABA_profile_center', GABA_export)


    n_e_it = 0
    n_i_it = 0
    AMPA_export = np.reshape(ei_net.AMPA_conn.W.todense()[n_e_it, :],
            (ei_net.Ni_y, ei_net.Ni_x))
    GABA_export = np.reshape(ei_net.GABA_conn1.W.todense()[n_i_it, :] +
        ei_net.extraGABA_conn1.W.todense()[n_i_it, :], (ei_net.Ne_y,
            ei_net.Ne_x))

    h5file.createArray(h5file.root, 'AMPA_profile_edges', AMPA_export)
    h5file.createArray(h5file.root, 'GABA_profile_edges', GABA_export)




    h5file.close()

