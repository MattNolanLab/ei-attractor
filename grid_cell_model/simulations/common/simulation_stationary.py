'''Main simulation run: Simulation of a stationary bump.'''
from __future__ import absolute_import, print_function, division

from numpy.random import choice
from nest.hl_api  import NESTError

from grid_cell_model.models.parameters  import getOptParser
from grid_cell_model.models.gc_net_nest import BasicGridCellNetwork
from grid_cell_model.data_storage       import DataStorage


parser          = getOptParser()
(options, args) = parser.parse_args()

output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
d = DataStorage.open(output_fname, 'a')
if ("trials" not in d.keys()):
    d['trials'] = []


overalT = 0.
################################################################################
for trial_idx in range(len(d['trials']), options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    d['invalidated'] = 1
    try:
        ei_net = BasicGridCellNetwork(options, simulationOpts=None)
        
        const_v = [0.0, 0.0]
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
        d['trials'].append(ei_net.getAllData())
        d.flush()
        constrT, simT, totalT = ei_net.printTimes()
        overalT += totalT
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
        break

d.close()
print("Script total run time: {0} s".format(overalT))
################################################################################

