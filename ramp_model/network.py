#
#   network.py
#
#   Network setup for Entorhinal cortical gamma constant/ramp input current
#   experiments
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

import numpy    as np
import logging  as lg
import nest

from scipy          import linspace
from numpy.random   import rand

from common         import *
from log            import *

__all__ = ['GammaNetwork']

class GammaNetwork:

    def initStates(self):
        # Initialize membrane potential randomly
        #self.E_pop.vm = self.EL_e + (self.Vt_e-self.EL_e) * rand(len(self.E_pop))
        #self.I_pop.vm = self.EL_i + (self.Vt_i-self.EL_i) * rand(len(self.I_pop))
        #self.setBackgroundInput(self.o.Iext_e*amp, self.o.Iext_i*amp)
        #self.E_pop.EL = (self.o.EL_e - self.o.EL_e_spread/2. +
        #        self.o.EL_e_spread*rand(len(self.E_pop))) * volt
        #self.E_pop.taum = (self.o.taum_e - self.o.taum_e_spread/2. +
        #        self.o.taum_e_spread*rand(len(self.E_pop))) * second
        #self.I_pop.tau_ad = (self.o.ad_tau_i_mean +
        #        self.o.ad_tau_i_std*np.random.randn(len(self.I_pop.tau_ad)))*second
        #self.I_pop.EL = (self.o.EL_i - self.o.EL_i_spread/2.0 +
        #        self.o.EL_i_spread*rand(len(self.I_pop))) * second
        #self.I_pop.taum = (self.o.taum_i - self.o.taum_i_spread/2. +
        #        self.o.taum_i_spread*rand(len(self.I_pop))) * second
        pass

    def reinit(self):
        self.net.reinit(states=True)
        self.initStates()

    def __init__(self, options):
        self.o = options
        self.externalDevicesSet = False

        nest.Install('gridcells')
        nest.ResetKernel()

        # Create the network using nest, do not connect anything yet
        nest.SetKernelStatus({"resolution": self.o.sim_dt, "print_time": self.o.print_time})
        nest.SetKernelStatus({"local_num_threads": self.o.numThreads})
        

        
        log_info('GammaNetwork', 'Building integrate and fire network using nest')
        C_m_e = self.o.taum_e * self.o.gL_e
        C_m_i = self.o.taum_i * self.o.gL_i
        
        neuron_params_e = {"V_m"              : self.o.EL_e,
                           "C_m"              : C_m_e,
                           "t_ref"            : self.o.t_ref_e,
                           "V_peak"           : self.o.V_peak_e,
                           "V_reset"          : self.o.Vr_e,
                           "E_L"              : self.o.EL_e,
                           "g_L"              : self.o.gL_e,
                           "Delta_T"          : self.o.deltaT_e,
                           "V_th"             : self.o.Vt_e,
                           "E_AMPA"           : self.o.E_AMPA,
                           "E_GABA_A"         : self.o.E_GABA_A,
                           "tau_AMPA_fall"    : self.o.tau_AMPA,
                           "tau_GABA_A_rise"  : self.o.tau_GABA_A_rise,
                           "tau_GABA_A_fall"  : self.o.tau_GABA_A_fall,
                           "I_e"              : 0.0}
        
        neuron_params_i = {"V_m"              : self.o.EL_i,
                           "C_m"              : C_m_i,
                           "t_ref"            : self.o.t_ref_i,
                           "V_peak"           : self.o.V_peak_i,
                           "V_reset"          : self.o.Vr_i,
                           "E_L"              : self.o.EL_i,
                           "g_L"              : self.o.gL_i,
                           "Delta_T"          : self.o.deltaT_i,
                           "V_th"             : self.o.Vt_i,
                           "E_AMPA"           : self.o.E_AMPA,
                           "E_GABA_A"         : self.o.E_GABA_A,
                           "tau_AMPA_fall"    : self.o.tau_AMPA,
                           "tau_GABA_A_rise"  : self.o.tau_GABA_A_rise,
                           "tau_GABA_A_fall"  : self.o.tau_GABA_A_fall,
                           "I_e"              : 0.0}
        
        
        self.model_name  = "iaf_gridcells"
        self.nodes_ex    = nest.Create(self.model_name, self.o.Ne, params = neuron_params_e)
        self.nodes_in    = nest.Create(self.model_name, self.o.Ni, params = neuron_params_i)
        self.receptors   = nest.GetDefaults(self.model_name)['receptor_types']

        self.initStates()


    def connRandom(self):
        '''
        Connect each pair of neurons with a random probability given by:
          
          E --> I given by AMPA_density <0, 1>
          I --> E given by GABA_density <0, 1>
        '''

        g_AMPA_mean = self.o.g_AMPA_total / self.o.Ne / self.o.AMPA_density
        g_AMPA_sigma = np.sqrt(np.log(1 + self.o.g_AMPA_std**2/g_AMPA_mean**2))
        g_AMPA_mu = np.log(g_AMPA_mean) - 1/2*g_AMPA_sigma**2
        g_GABA_A_mean = self.o.g_GABA_total / self.o.Ni / self.o.GABA_density

        #self.AMPA_conn.connect_random(self.E_pop, self.I_pop, AMPA_density, 
        #        weight=lambda: np.random.lognormal(g_AMPA_mu, g_AMPA_sigma)*siemens)
        #self.GABA_conn1.connect_random(self.I_pop, self.E_pop, GABA_density,
        #        weight=self.B_GABA*g_GABA_mean)
        #self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W) 

        nest.CopyModel("static_synapse", "ex_AMPA",
                {"weight"        : g_AMPA_mean,
                 "delay"         : self.o.delay,
                 "receptor_type" : self.receptors["AMPA"]})
        nest.CopyModel("static_synapse", "in_GABA_A",
                {"weight"        : g_GABA_A_mean,
                 "delay"         : self.o.delay,
                 "receptor_type" : self.receptors["GABA_A"]})

        CE = self.o.Ne * self.o.AMPA_density
        CI = self.o.Ni * self.o.GABA_density
        nest.RandomConvergentConnect(self.nodes_ex, self.nodes_in, int(CE), model="ex_AMPA")
        nest.RandomConvergentConnect(self.nodes_in, self.nodes_ex, int(CI), model="in_GABA_A")
  


    ############################################################################ 
    #                     External sources definitions
    ############################################################################ 
    def setRampCurrent(self):
        '''
        Set up the ramp current in the model. Appropriate parameters from
        neuronOpts will be used here:
            Iext_start
            Iext_e_max
            Iext_i_max
        '''
        nest.SetStatus(self.nodes_ex, {"I_e" : 0.0})
        nest.SetStatus(self.nodes_in, {"I_e" : 0.0})
        model_name = 'ramp_current_generator'

        # E ramp current
        params_e = {
                'a' :  self.o.Iext_e_max / (self.o.time - self.o.Iext_start),
                'b' : -self.o.Iext_e_max / (self.o.time - self.o.Iext_start) * self.o.Iext_start,
                'c' : 1.0}
        params_i = {
                'a' :  self.o.Iext_i_max / (self.o.time - self.o.Iext_start),
                'b' : -self.o.Iext_i_max / (self.o.time - self.o.Iext_start) * self.o.Iext_start,
                'c' : 1.0}
        self._nodes_ramp_e = nest.Create(model_name, self.o.Ne, params=params_e)
        self._nodes_ramp_i = nest.Create(model_name, self.o.Ni, params=params_i)

        nest.Connect(self._nodes_ramp_e, self.nodes_ex, model="static_synapse")
        nest.Connect(self._nodes_ramp_i, self.nodes_in, model="static_synapse")

        log_warn('GammaNetwork', 'Ramp current inputs need the Gaussian profile')

    def setConstantCurrent(self):
        '''
        Set up constant external current excitation, given by parameters in
        neuronOpts:
            Iext_e
            Iext_i
        '''
        if self.externalDevicesSet == True:
            log_warn('ramp_model.GammaNetwork', 'Setting external excitation sources twice might be a bug!')
        else:
            self.externalDevicesSet = True

        nest.SetStatus(self.nodes_ex, {"I_e" : self.o.Iext_e})
        nest.SetStatus(self.nodes_in, {"I_e" : self.o.Iext_i})

    def setRampPoissonConductance(self):
        '''
        Set up external conductance connections, which are Poisson processes.
        Each neuron will receive an independent poisson ramp conductance, given
        by parameters from neuronOpts:
            g_ext_e
            g_ext_i
            rate_ext_e_min
            rate_ext_e_max
            rate_ext_i_min
            rate_ext_i_max
        '''
        raise NotImplementedException("GammaNetwork.setRampPoissonConductance")

    def setConstantPoissonConductance(self):
        '''
        Set a constant external poisson conductance, given by these parameters
        in neuronOpts:
            g_ext_e
            g_ext_i
            rate_ext_e
            rate_ext_i
        ''' 
        raise NotImplementedException("GammaNetwork.setConstantPoissonConductance")


    ############################################################################ 
    #                              Monitors
    ############################################################################ 
    def monitorESpikes(self, nids):
        '''
        Monitor excitatory neurons' spikes:
          nids  - neuron ids (local, within the group)
        '''
        self._espikes = nest.Create("spike_detector")
        self._sp_e_nids = nids
        nest.SetStatus([self._espikes],[{"label": "spikes_ex",
                           "withtime": True,
                           "withgid": True}])
        nest.ConvergentConnect(list(self.nodes_ex[0] + np.array(nids)), self._espikes, model="static_synapse")

    def monitorISpikes(self, nids):
        '''
        Monitor inhibitory neurons' spikes:
          nids  - neuron ids (local, within the group)
        '''
        self._ispikes = nest.Create("spike_detector")
        self._sp_i_nids = nids
        nest.SetStatus([self._ispikes],[{"label": "ramp_model_in",
                   "withtime": True,
                   "withgid": True}])
        nest.ConvergentConnect(list(self.nodes_in[0] + np.array(nids)), self._ispikes, model="static_synapse")


    def getESpikeMonitor(self):
        return self._espikes, self._sp_e_nids


    def getISpikeMonitor(self):
        return self._ispikes, self._sp_i_nids


    def monitorEState(self, nids, states):
        '''
        Monitor state variables of neurons in the E population.
          nids  - list of ids of neurons within the E population (starting at 0)
          state - a list of state names to record
        '''
        self._state_e_nids  = nids
        self._states_e      = states
        self._meter_e       = nest.Create('multimeter',
                                params = {'withtime': True,
                                          'interval': 0.1,
                                          'record_from': states})
        nest.Connect(self._meter_e, np.array(nids) + self.nodes_ex[0])

    def monitorIState(self, nids, states):
        '''
        Monitor state variables of neurons in the I population.
          nids  - list of ids of neurons within the E population (starting at 0)
          state - a list of state names to record
        '''
        self._state_i_nids  = nids
        self._states_i      = states
        self._meter_i       = nest.Create('multimeter',
                                params = {'withtime': True,
                                          'interval': 0.1,
                                          'record_from': states})
        nest.Connect(self._meter_i, np.array(nids) + self.nodes_in[0])


