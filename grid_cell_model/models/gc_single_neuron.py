#
#   gc_single_neuron.py
#
#   A model that simulates a single neuron with theta input. One neuron from E
#   and I population.
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

import numpy    as np
import logging  as lg

from numpy.random import rand, randn

import gc_neurons
from gc_net       import GridCellNetwork
from data_storage import DataStorage

import nest
nest.Install('gridcellsmodule')

class OneNeuronNetwork(GridCellNetwork):
    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        # Artificially change net_Ne and net_Ni so that we can use only 1
        # neuron per each group
        self.net_Ne = 1
        self.net_Ni = 1
        self.Ne_x = 1
        self.Ne_y = 1
        self.Ni_x = 1
        self.Ni_y = 1


        self.spikeMon_e = None
        self.spikeMon_i = None
        self.stateMon_e = None
        self.stateMon_i = None

        self._initNESTKernel()
        self._constructNetwork()
        self._initStates()
        self._initCellularProperties()

        # State monitors
        self.stateMonParams_e = self.getDefaultStateMonParams()
        self.stateMonParams_i = self.getDefaultStateMonParams()

        self.stateMon_e  = self.getStateMonitor("E", [0], self.stateMonParams_e)
        self.stateMon_i  = self.getStateMonitor("I", [0], self.stateMonParams_i)


    def getDefaultStateMonParams(self):
        return {
            'withtime' : True,
            'interval' : self.no.sim_dt,
            'record_from' : ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
                'I_clamp_GABA_A', 'I_stim']
        }



    def uniformDistrib(self, mean, spread, N):
        return mean - spread/2.0 * rand(N)

    def _initStates(self):
        # Initialize membrane potential randomly
        nest.SetStatus(self.E_pop, 'V_m', self.no.EL_e + (self.no.Vt_e-self.no.EL_e) * np.random.rand(len(self.E_pop)))
        nest.SetStatus(self.I_pop, 'V_m', self.no.EL_i + (self.no.Vt_i-self.no.EL_i) * np.random.rand(len(self.I_pop)))


    def _initCellularProperties(self):
        EL_e    = self.uniformDistrib(self.no.EL_e,   self.no.EL_e_spread,   len(self.E_pop))
        taum_e  = self.uniformDistrib(self.no.taum_e, self.no.taum_e_spread, len(self.E_pop))
        EL_i    = self.uniformDistrib(self.no.EL_i,   self.no.EL_i_spread,   len(self.I_pop))
        taum_i  = self.uniformDistrib(self.no.taum_i, self.no.taum_i_spread, len(self.I_pop))
        nest.SetStatus(self.E_pop, 'E_L',  EL_e)
        nest.SetStatus(self.E_pop, 'C_m',  taum_e * self.no.gL_e)
        nest.SetStatus(self.I_pop, 'E_L',  EL_i)
        nest.SetStatus(self.I_pop, 'C_m', taum_i * self.no.gL_i)

        #self.I_pop.tau_ad   = (self.no.ad_tau_i_mean + self.no.ad_tau_i_std * np.random.randn(len(self.I_pop.tau_ad))) * ms


    def _initNESTKernel(self):
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution" : self.no.sim_dt, "print_time": False})
        nest.SetKernelStatus({"local_num_threads" : self.no.nthreads})


    def _constructNetwork(self):
        '''Construct the E/I network'''
        self.e_neuron_params = gc_neurons.getENeuronParams(self.no)
        self.i_neuron_params = gc_neurons.getINeuronParams(self.no)


        self.e_model_name = "iaf_gridcells"
        self.i_model_name = "iaf_gridcells"
        self.E_pop = nest.Create(self.e_model_name, self.net_Ne,
                params=self.e_neuron_params)
        self.I_pop = nest.Create(self.i_model_name, self.net_Ni, params =
                self.i_neuron_params)



    def simulate(self, time, printTime=True):
        '''Run the simulation'''
        nest.SetKernelStatus({"print_time": bool(printTime)})
        nest.Simulate(time)


    def getStateMonitor(self, type, N_ids, params):
        '''
        Return a state monitor for a given population (type) and relative
        indexes of neurons (N_ids), with parameters given by params
        '''
        if (len(N_ids) == 0):
            raise ValueError("State monitor needs to record from at least one neuron")

        N = len(N_ids)

        if (type == "E"):
            if self.stateMon_e is None:
                self.stateMon_e = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_e, self.E_pop[0] + np.array(N_ids))
            return self.stateMon_e
        elif (type == "I"):
            if (self.stateMon_i is None):
                self.stateMon_i = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_i, self.I_pop[0] + np.array(N_ids))
            return self.stateMon_i


    def getAttrDictionary(self):
        d = {}

        d['e_neuron_params'] = self.e_neuron_params
        d['i_neuron_params'] = self.i_neuron_params
        d['E_pop'          ] = np.array(self.E_pop)
        d['I_pop'          ] = np.array(self.I_pop)
        d['net_Ne'         ] = self.net_Ne
        d['net_Ni'         ] = self.net_Ni
        d['Ne_x'           ] = self.Ne_x
        d['Ne_y'           ] = self.Ne_y
        d['Ni_x'           ] = self.Ni_x
        d['Ni_y'           ] = self.Ni_y

        return d
            

    def getNetParams(self):
        out = {}
        out['options']  = self.no._einet_optdict
        out['net_attr'] = self.getAttrDictionary()
        return out


    def saveData(self):
        output_fname = "{0}/{1}noise_sigma{2}_output.h5".format(self.no.output_dir,
                self.no.fileNamePrefix, int(self.no.noise_sigma))

        out = self.getNetParams()
        out['stateMon_e']   = nest.GetStatus(self.stateMon_e)
        out['stateMon_i']   = nest.GetStatus(self.stateMon_i)

        d = DataStorage.open(output_fname, 'w')
        d.update(out)
        d.close()




