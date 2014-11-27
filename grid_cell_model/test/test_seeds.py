'''Integration tests of reproducibility with random seeds.'''
from __future__ import absolute_import, print_function, division

import pytest
import numpy as np

from nest.hl_api        import NESTError
from grid_cell_model.models.gc_net_nest import BasicGridCellNetwork, ConstPosInputs
from grid_cell_model.models.seeds import TrialSeedGenerator

from data.network_params import defaultParameters

class Params(object):
    '''Dictionary that has fields as natural naming.'''
    def __init__(self, params):
        self._params = params
        self.__dict__ = params

def run_network(options):
    '''Run network with given `options` and return the dictionary that would be
    saved to an output file.
    '''
    o = options

    seed_gen = TrialSeedGenerator(o.master_seed)

    stateMonParams = {
            'start' : o.time - o.stateMonDur
    }
    nrec_spikes_e = None # all neurons
    nrec_spikes_i = 10

    d = {}
    if "trials" not in d.keys():
        d['trials'] = []

    overalT = 0.
    ################################################################################
    for trial_idx in range(len(d['trials']), o.ntrials):
        print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
        seed_gen.set_generators(trial_idx)
        d['invalidated'] = 1
        ei_net = BasicGridCellNetwork(o, simulationOpts=None,
                nrec_spikes=(nrec_spikes_e, nrec_spikes_i),
                stateRecParams=(stateMonParams, stateMonParams))
        if o.velON and not o.constantPosition:
            ei_net.setVelocityCurrentInput_e()
        posIn = ConstPosInputs(0, 0) if o.constantPosition else None
        ei_net.setPlaceCells(posIn=posIn)

        try:
            ei_net.simulate(o.time, printTime=o.printTime)
        except NESTError as e:
            print("Simulation interrupted. Message: {0}".format(str(e)))
            print("Trying to save the simulated data if possible...")
        ei_net.endSimulation()
        d['trials'].append(ei_net.getAllData())
        constrT, simT, totalT = ei_net.printTimes()
        overalT += totalT

    print("Script total run time: {0} s".format(overalT))
    ################################################################################
    return d


def cmp_dict(d1, d2):
    '''A simple comparator that assumes number of items are the same.'''
    if len(d1) != len(d2):
        return False
    for key in d1.keys():
        print("key: %s" % key)

        if not key in d2.keys():
            return False

        if isinstance(d1[key], dict):
            if not cmp_dict(d1[key], d2[key]):
                return False
        elif isinstance(d1[key], list):
            if not cmp_list(d1[key], d2[key]):
                return False
        elif isinstance(d1[key], np.ndarray):
            if not np.all(d1[key] == d2[key]):
                return False
        else:
            if d1[key] != d2[key]:
                return False
    return True


def cmp_list(l1, l2):
    '''A simple comparator that compares two lists recursively.'''
    if len(l1) != len(l2):
        return False
    i = 0
    for val1, val2 in zip(l1, l2):
        print("i: %d" % i)
        if isinstance(val1, dict):
            if not cmp_dict(val1, val2):
                return False
        elif isinstance(val1, list):
            if not cmp_list(val1, val2):
                return False
        elif isinstance(val1, np.ndarray):
            if not np.all(val1 == val2):
                return False
        else:
            if val1 != val2:
                return False

        i += 1
    return True


@pytest.fixture(scope='function')
def fix_params():
    return defaultParameters


def test_short_animal_movement(fix_params):
    '''Test several trials of short simulations of animal movement.'''
    p = fix_params.copy()
    p['noise_sigma'] = 150.0 # pA

    # Submitting
    p['time']             = 1e3  # ms
    p['nthreads']         = 1
    p['ntrials']          = 2
    p['velON']            = 1
    p['constantPosition'] = 0
    p['master_seed']      = 123456

    p['_einet_optdict'] = p.copy()

    data1 = run_network(Params(p))
    data2 = run_network(Params(p))
    #test_data = {
    #    "one" : [
    #        {"one1": 10, "one2": np.arange(10), "one3": [0, 1, 2]},
    #        {"one1": 11, "one2": np.arange(100)},
    #    ],
    #}
    assert cmp_dict(data1, data2) is True
