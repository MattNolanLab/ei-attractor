from brian import *
from brian.library.IF import *
from brian.library.synapses import *
from brian.membrane_equations import *

from scipy import linspace
from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

import numpy as np

import time
import math
import random



class EI_Network:
    def __init__(o, clk):
    # Setup neuron equations
    # Using exponential integrate and fire model
    # Excitatory population
    Rm_e = o.Rm_e*Mohm
    taum_e = o.taum_e*msecond
    Ce = taum_e/Rm_e
    EL_e = o.EL_e*mV
    deltaT_e = o.deltaT_e
    Vt_e = o.Vt_e*mV
    Vr_e = o.Vr_e*mV
    tau_GABA_rise = o.tau_GABA_rise*msecond
    tau_GABA_fall = o.tau_GABA_fall*msecond
    tau1_GABA = tau_GABA_fall
    tau2_GABA = tau_GABA_rise*tau_GABA_fall / (tau_GABA_rise*tau_GABA_fall);
    B_GABA = 1/((tau2_GABA/tau1_GABA)^(tau_GABA_rise/tau1_GABA) - 
            (tau2_GABA/tau1_GABA)^(tau_GABA_rise/tau2_GABA))
    tau_ad_e = o.ad_tau_e_mean*msecond
    
    Vrev_GABA = o.Vrev_GABA*mV


    self.eqs_e = Equation('''
        dvm/dt = 1/C*Im + (noise_sigma*xi/taum**.5): volt
        Im = gL*(EL-vm)*(1+g_ad/gL)+gL*deltaT*exp((vm-Vt)/DeltaT) + (gi1 - gi2)*(Esyn - vm) + Iext  : amp
        dgi1/dt = -gi1/syn_tau1 : siemens
        dgi2/dt = -gi2/syn_tau2 : siemens
        dg_ad/dt = -g_ad/tau_ad : siemens
        Iext : amp
        ''',
        C=Ce,
        gL=1/Rm_e,
        EL=EL_e,
        deltaT=deltaT_e,
        Vt=Vt_e,
        Esyn=Vrev_GABA,
        syn_tau1=tau1_GABA,
        syn_tau2=tau2_GABA,
        tau_ad=tau_ad_e)


    # Inhibitory population
    Rm_i = o.Rm_i*Mohm
    taum_i = o.taum_i
    Ci = taum_i/Rm_i
    EL_i = o.EL_i.mV
    deltaT_i = o.deltaT_i*mV
    Vt_i = o.Vt_i*mV
    Vr_i = o.Vr_i*mV
    Vrev_AMPA = o.Vrev_AMPA*mV
    tau_AMPA = o.tau_AMPA*msecond
    tau_ad_i = o.ad_taui_mean*msecond
    
    self.eqs_i = Equation('''
        dvm/dt = 1/C*Im + (noise_sigma*xi/taum**.5): volt
        Im = gL*(EL-vm)*(1+g_ad/gL)+gL*deltaT*exp((vm-Vt)/DeltaT) + ge*(Esyn - vm) + Iext  : amp
        dge/dt = -ge/syn_tau : siemens
        dg_ad/dt = -g_ad/tau_ad : siemens
        Iext : amp
        ''',
        C=Ci,
        gL=1/Rm_i,
        EL=EL_i,
        deltaT=deltaT_i,
        Vt=Vt_i,
        Esyn=Vrev_AMPA,
        syn_tau=tau_AMPA,
        tau_ad=tau_ad_i)


    # Other constants
    refrac_abs = o.refrac_abs*msecond
    spike_detect_th = o.spike_detect_th*mV

    g_AMPA_sigma = np.sqrt(np.log(1 + o.g_AMPA_std**2/o.g_AMPA_mean**2))
    g_AMPA_mu = np.log(o.g_AMPA_mean) - 1/2*g_AMPA_sigma**2



    # Setup neuron groups and connections
    self.E_pop = NeuronGroup(
            N = o.Ne
            model=eqs_e,
            threshold=spike_detect_th,
            reset=Vr_e,
            refractory=refrac_abs,
            clock=clock)

    self.I_pop = NeuronGroup(
            N = o.Ni,
            model=eqs_i,
            threshold=spike_detect_th,
            reset=Vr_i,
            refractory=refrac_abs,
            clock=clock)

    self.AMPA_conn = Connection(E_pop, I_pop, 'ge')
    self.AMPA_conn.connect_random(E_pop, I_pop, o.AMPA_density,
            weight=lambda:np.random.lognormal(g_AMPA_mu, g_AMPA_sigma)*nS)

    self.GABA_conn1 = Connection(I_pop, E_pop, 'gi1')
    self.GABA_conn1.connect_random(I_pop, E_pop, o.GABA_density,
            weight=o.g_GABA_mean*nS))
    self.GABA_conn2 = Connection(I_pop, E_pop, 'gi2')
    self.GABA_conn2.connect(I_pop, E_pop, GABA_conn1.W)



    # Setup adaptation connections: neuron on itself
    self.adaptConn_e = IdentityConnection(E_pop, E_pop,  'g_ad',
            weight=o.ad_e_g_inc*siemens)
    self.adaptConn_i = IdentityConnection(I_pop, I_pop, 'g_ad',
            weight=o.ad_i_g_inc*siemens)



    # Initialize membrane potential randomly
    if (noise_sigma == 0):
        self.E_pop.vm = EL_e + (Vt_e-EL_e) * rand(len(self.E_pop))
        self.I_pop.vm = EL_i + (Vt_i-EL_i) * rand(len(self.I_pop))
    else:
        self.E_pop.vm = EL_e + zeros(len(self.E_pop))
        self.I_pop.vm = EL_i + zeros(len(self.I_pop))


    # Create network (withough any monitors
    self.net = Network(
            self.E_pop,
            self.I_pop,
            self.AMPA_conn,
            self.GABA_conn1,
            self.GABA_conn2,
            self.adaptConn_e,
            self.adaptConn_i)

