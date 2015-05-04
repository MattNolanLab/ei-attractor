'''Main simulation run: Simulation of a stationary bump.'''
from __future__ import absolute_import, print_function, division

from numpy.random import choice
from nest.hl_api import NESTError

from grid_cell_model.models.parameters import getOptParser
from grid_cell_model.models.gc_net_nest import BasicGridCellNetwork
from grid_cell_model.models.seeds import TrialSeedGenerator
from grid_cell_model.parameters.data_sets import DictDataSet
from grid_cell_model.visitors.spikes import SpikeStatsVisitor
from grid_cell_model.visitors.signals import AutoCorrelationVisitor
from simtools.storage import DataStorage


def signal_analysis(data):
    '''Run the signal analysis visitors on a single data trial.

    Parameters
    ----------
    data : dict
        A dictionary containing data of one trial.

    Returns
    -------
    data : dict
        Input data modified in-situ.
    '''
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']
    dummy_data_set = DictDataSet(data)

    stats_visitor_e = SpikeStatsVisitor("spikeMon_e", forceUpdate=False)
    ac_visitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=False)
    stats_visitor_e.visitDictDataSet(dummy_data_set)
    ac_visitor.visitDictDataSet(dummy_data_set)

    # Clean the state monitor
    data['stateMonF_e'] = [data['stateMonF_e'][0]]

    return data


parser          = getOptParser()
(options, args) = parser.parse_args()

output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
                                                   options.fileNamePrefix,
                                                   options.job_num)
d = DataStorage.open(output_fname, 'a')
if "trials" not in d.keys():
    d['trials'] = []

seed_gen = TrialSeedGenerator(options.master_seed)

overalT = 0.

###############################################################################
for trial_idx in range(len(d['trials']), options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    seed_gen.set_generators(trial_idx)
    d['master_seed'] = options.master_seed
    d['invalidated'] = 1
    try:
        ei_net = BasicGridCellNetwork(options, simulationOpts=None)

        const_v = [0.0, 0.0]
        ei_net.setConstantVelocityCurrent_e(const_v)

        stateRecF_e = choice(ei_net.E_pop, options.gammaNSample, replace=False)
        stateRecF_i = choice(ei_net.I_pop, options.gammaNSample, replace=False)

        stateMonF_e_params = {
            'withtime': False,
            'interval': options.sim_dt * 10,
            'record_from': ['I_clamp_GABA_A']
        }

        stateMonF_e = ei_net.getGenericStateMonitor(stateRecF_e,
                                                    stateMonF_e_params,
                                                    'stateMonF_e')

        d['net_params'] = ei_net.getNetParams()  # Common settings will stay
        d.flush()

        ei_net.simulate(options.time, printTime=options.printTime)
        ei_net.endSimulation()
        d['trials'].append(signal_analysis(ei_net.getAllData()))
        d.flush()
        constrT, simT, totalT = ei_net.printTimes()
        overalT += totalT
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
        break

d.close()
print("Script total run time: {0} s".format(overalT))
###############################################################################
