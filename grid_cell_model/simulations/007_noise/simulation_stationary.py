#
#   simulation_stationary.py
#
#   Main simulation run: Simulation of a stationary bump.
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
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork
from data_storage       import DataStorage


parser          = getOptParser()
(options, args) = parser.parse_args()


out = []
overalT = 0.
################################################################################
for trial_idx in range(options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    ei_net = BasicGridCellNetwork(options, simulationOpts=None)
    
    const_v = [00.0, 0.0]
    ei_net.setConstantVelocityCurrent_e(const_v)
    
    
    stateRecF_e = choice(ei_net.E_pop, options.gammaNSample, replace=False)
    stateRecF_i = choice(ei_net.I_pop, options.gammaNSample, replace=False)
    
    stateMonF_e_params = {
            'withtime' : False,
            'interval' : options.sim_dt*10,
            'record_from' : ['I_clamp_GABA_A']
    }
    stateMonF_i_params = dict(stateMonF_e_params)
    stateMonF_i_params['record_from'] = ['I_clamp_AMPA', 'I_clamp_NMDA']

    stateMonF_e = ei_net.getGenericStateMonitor(stateRecF_e,
            stateMonF_e_params, 'stateMonF_e')
    stateMonF_i = ei_net.getGenericStateMonitor(stateRecF_i,
            stateMonF_i_params, 'stateMonF_i')
    
    ei_net.simulate(options.time, printTime=options.printTime)
    ei_net.endSimulation()
    out.append(ei_net.getAllData())
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT

output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
d = DataStorage.open(output_fname, 'w')
d["trials"] = out
d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################

