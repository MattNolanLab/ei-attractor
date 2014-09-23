'''
.. currentmodule: test.test_single_neuron

The :mod:`~test.test_single_neuron` module tests the basic functionality of
single neurons.
'''
from __future__ import absolute_import, division, print_function

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pytest
import numpy as np
import nest

from grid_cell_model.models.gc_single_neuron import NeuronAndGenerator
from data import single_neuron_params

@pytest.fixture(scope='module', params=[.0, 0.1, 1.0])
def C_Mg_fix(request):
    return request.param


@pytest.fixture(scope='module')
def default_params_fix(C_Mg_fix):
    class NeuronOptions(object):
        def __init__(self, d):
            self.no = d
            self.__dict__ = d

    p = single_neuron_params.defaultParameters.copy()
    p['C_Mg'] = C_Mg_fix
    return NeuronOptions(p)


@pytest.fixture(scope='module')
def T_fix():
    return 1e3


def nmda_gating_testgen(Vm, mg_c):
    '''cf. Jahr & Stevens (1990)'''
    return 1. / (1 + mg_c / 3.57 * np.exp(-0.062 * Vm) )


class TestNMDA(object):

    def _simulate(self, params, T):
        network = NeuronAndGenerator(params, simulationOpts=None)
        network.simulate(T, printTime=False)
        return network

    def _extract_state_vars(self, state_mon):
        stat = nest.GetStatus(state_mon)
        return stat[0]['events']

    def test_nmda_gating_variable(self, default_params_fix, T_fix):
        network = self._simulate(default_params_fix, T_fix)
        vars = self._extract_state_vars(network.stateMon_i)

        np.testing.assert_allclose(
                vars['s_NMDA'],
                nmda_gating_testgen(vars['V_m'], network.no.C_Mg),
                rtol=1e-10)

    def test_examine_state_variables(self, default_params_fix, T_fix):
        default_params_fix.Iext_i_theta = 0.

        network = self._simulate(default_params_fix, T_fix)
        e = self._extract_state_vars(network.stateMon_i)

        N = 4
        fig = plt.figure()
        ax1 = fig.add_subplot(N, 1, 1)
        ax1.plot(e['times'], e['V_m'])
        ax1.set_ylabel('Vm (mV)')

        ax2 = fig.add_subplot(N, 1, 2)
        ax2.plot(e['times'], e['s_NMDA'])
        ax2.set_ylabel('s_NMDA')
        ax2.set_xlabel('t (ms)')

        ax3 = fig.add_subplot(N, 1, 3)
        ax3.plot(e['times'], e['g_NMDA'])
        ax3.set_ylabel('g_NMDA (nS)')
        ax3.set_xlabel('t (ms)')

        ax4 = fig.add_subplot(N, 1, 4)
        ax4.plot(e['times'], e['I_stim'])
        ax4.set_ylabel('I_stim (pA)')
        ax4.set_xlabel('t (ms)')

        fname = 'single_neuron_state_variables_Mg%.1f_mM.pdf' % network.no.C_Mg
        fig.savefig(fname)


