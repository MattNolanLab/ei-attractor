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

#def AMPA_p(i, j, mu, sigma):
#    return np.exp(-(abs(i-j) - mu)**2/2.0/sigma**2)
#
#def GABA_p(i, j, sigma):
#    return np.exp(-(abs(i-j))**2/2.0/sigma**2)


class EI_Network:
    def __init__(self, o, clk):

        noise_sigma=o.noise_sigma*volt

        # Setup neuron equations
        # Using exponential integrate and fire model
        # Excitatory population
        Rm_e = o.Rm_e*ohm
        taum_e = o.taum_e*second
        Ce = taum_e/Rm_e
        EL_e = o.EL_e*volt
        deltaT_e = o.deltaT_e*volt
        Vt_e = o.Vt_e*volt
        Vr_e = o.Vr_e*volt
        tau_GABA_rise = o.tau_GABA_rise*second
        tau_GABA_fall = o.tau_GABA_fall*second
        tau1_GABA = tau_GABA_fall
        tau2_GABA = tau_GABA_rise*tau_GABA_fall / (tau_GABA_rise + tau_GABA_fall);
        B_GABA = 1/((tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau1_GABA) - 
                (tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau2_GABA))
        tau_ad_e = o.ad_tau_e_mean*second
        
        Vrev_GABA = o.Vrev_GABA*volt


        self.eqs_e = Equations('''
            dvm/dt = 1/C*Im + (noise_sigma*xi/taum**.5): volt
            Im = gL*(EL-vm)*(1+g_ad/gL)+gL*deltaT*exp((vm-Vt)/deltaT) + Isyn + Iext  : amp
            Isyn = (gi1 - gi2)*(Esyn - vm) : amp
            dgi1/dt = -gi1/syn_tau1 : siemens
            dgi2/dt = -gi2/syn_tau2 : siemens
            dg_ad/dt = -g_ad/tau_ad : siemens
            Iext : amp
            ''',
            C=Ce,
            gL=1/Rm_e,
            noise_sigma=noise_sigma,
            taum=taum_e,
            EL=EL_e,
            deltaT=deltaT_e,
            Vt=Vt_e,
            Esyn=Vrev_GABA,
            syn_tau1=tau1_GABA,
            syn_tau2=tau2_GABA,
            tau_ad=tau_ad_e)


        # Inhibitory population
        Rm_i = o.Rm_i*ohm
        taum_i = o.taum_i*second
        Ci = taum_i/Rm_i
        EL_i = o.EL_i*volt
        deltaT_i = o.deltaT_i*volt
        Vt_i = o.Vt_i*volt
        Vr_i = o.Vr_i*volt
        Vrev_AMPA = o.Vrev_AMPA*volt
        tau_AMPA = o.tau_AMPA*second
        tau_ad_i = o.ad_tau_i_mean*second
        
        self.eqs_i = Equations('''
            dvm/dt = 1/C*Im + (noise_sigma*xi/taum**.5): volt
            Im = gL*(EL-vm)*(1+g_ad/gL) + gL*deltaT*exp((vm-Vt)/deltaT) + Isyn + Iext  : amp
            Isyn = ge*(Esyn - vm) : amp
            dge/dt = -ge/syn_tau : siemens
            dg_ad/dt = -g_ad/tau_ad : siemens
            Iext : amp
            ''',
            C=Ci,
            gL=1/Rm_i,
            noise_sigma=noise_sigma,
            taum=taum_i,
            EL=EL_i,
            deltaT=deltaT_i,
            Vt=Vt_i,
            Esyn=Vrev_AMPA,
            syn_tau=tau_AMPA,
            tau_ad=tau_ad_i)


        # Other constants
        refrac_abs = o.refrac_abs*second
        spike_detect_th = o.spike_detect_th*volt

        g_AMPA_mean = o.g_AMPA_total/o.Ne
        g_AMPA_sigma = np.sqrt(np.log(1 + o.g_AMPA_std**2/g_AMPA_mean**2))
        g_AMPA_mu = np.log(g_AMPA_mean) - 1/2*g_AMPA_sigma**2
        #g_AMPA_mean = o.g_AMPA_mean * siemens
        g_GABA_mean = o.g_GABA_total / o.Ni * siemens



        # Setup neuron groups and connections
        self.E_pop = NeuronGroup(
                N = o.Ne,
                model=self.eqs_e,
                threshold=spike_detect_th,
                reset=Vr_e,
                refractory=refrac_abs,
                clock=clk)

        self.I_pop = NeuronGroup(
                N = o.Ni,
                model=self.eqs_i,
                threshold=spike_detect_th,
                reset=Vr_i,
                refractory=refrac_abs,
                clock=clk)

        # Generate connection-probability profile functions for GABA and AMPA connections
        self.pAMPA_mu = 0.4 * o.Ni
        self.pAMPA_sigma = 0.5/6 * o.Ni
        self.pGABA_sigma = 0.25/6 * o.Ne

        self.generate_pAMPA_template(o.Ni, self.pAMPA_mu, self.pAMPA_sigma)
        self.generate_pGABA_template(o.Ne, self.pGABA_sigma)

        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge')
        for i in xrange(o.Ne):
            e_norm = int(round(np.double(i)/o.Ne*o.Ni))

            tmp_templ = np.roll(self.pAMPA_templ, e_norm)
            self.AMPA_conn.W.rows[i] = list((rand(o.Ni) <
                o.AMPA_density*tmp_templ).nonzero()[0])
            self.AMPA_conn.W.data[i] = np.random.lognormal(g_AMPA_mu,
                    g_AMPA_sigma, len(self.AMPA_conn.W.rows[i]))*siemens

        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        for i in xrange(o.Ni):
            i_norm = int(round(np.double(i)/o.Ni*o.Ne))
            

            tmp_templ = np.roll(self.pGABA_templ, i_norm)
            self.GABA_conn1.W.rows[i] = list((rand(o.Ne) <
                o.GABA_density*tmp_templ).nonzero()[0])
            self.GABA_conn1.W.data[i] = [B_GABA*g_GABA_mean] * len(self.GABA_conn1.W.rows[i])

        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)



        # Setup adaptation connections: neuron on itself
        self.adaptConn_e = IdentityConnection(self.E_pop, self.E_pop,  'g_ad',
                weight=o.ad_e_g_inc*siemens)
        self.adaptConn_i = IdentityConnection(self.I_pop, self.I_pop, 'g_ad',
                weight=o.ad_i_g_inc*siemens)



        # Initialize membrane potential randomly
        self.E_pop.vm = EL_e + (Vt_e-EL_e) * rand(len(self.E_pop))
        self.I_pop.vm = EL_i + (Vt_i-EL_i) * rand(len(self.I_pop))


        # Create network (withough any monitors
        self.net = Network(
                self.E_pop,
                self.I_pop,
                self.AMPA_conn,
                self.GABA_conn1,
                self.GABA_conn2,
                self.adaptConn_e,
                self.adaptConn_i)

        self.setBackgroundInput(o.Iext_e*amp, o.Iext_i*amp)

        self.o = o


    def setBackgroundInput(self, Iext_e, Iext_i):
        self.E_pop.Iext = linspace(Iext_e, Iext_e, len(self.E_pop))
        self.I_pop.Iext = linspace(Iext_i, Iext_i, len(self.I_pop))

    def generate_pAMPA_template(self, N, mu, sigma):
        '''Generate AMPA probability profile function on an interval [0, N-1]. For
        now it will be a Gaussian profile, with mean mu and std. dev. sigma. The
        boundaries are wrapped-around (N==0). This is distance-dependent
        profile, which take wrap around boundaries into account.'''
        self.pAMPA_templ = np.exp(-(np.arange(N/2+1) - mu)**2/2.0/sigma**2)
        self.pAMPA_templ = np.concatenate((self.pAMPA_templ,
            self.pAMPA_templ[1:N-len(self.pAMPA_templ)+1][::-1]))
        
    def generate_pGABA_template(self, N, sigma):
        '''Generate GABA probability profile function on an interval [0, N-1].
        It is similar to generate_pAMPA_template but the mean is 0'''
        self.pGABA_templ = np.exp(-(np.arange(N/2+1))**2/2.0/sigma**2)
        self.pGABA_templ = np.concatenate((self.pGABA_templ,
            self.pGABA_templ[1:N-len(self.pGABA_templ)+1][::-1]))


