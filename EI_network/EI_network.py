from brian import *
from brian.library.IF import *
from brian.library.synapses import *
from brian.membrane_equations import *

from scipy import linspace
from scipy.io import loadmat
from optparse import OptionParser
from datetime import datetime

import numpy as np
from numpy.random import rand

import time
import math



class EI_Network:
    def __init__(self, o, clk):
        if o.ndim == 1 or o.ndim == -1:
            self.net_Ne = o.Ne
            self.net_Ni = o.Ni
        elif o.ndim == 2:
            self.net_Ne = o.Ne**2
            self.net_Ni = o.Ni**2
        else:
            raise Exception("Number of Mexican hat dimensions must be 0 or 1" +
                ", not" + str(ndim) + ".")

        self._createBasicNet(o, clk)

    def initStates(self):
        # Initialize membrane potential randomly
        self.E_pop.vm = self.EL_e + (self.Vt_e-self.EL_e) * rand(len(self.E_pop))
        self.I_pop.vm = self.EL_i + (self.Vt_i-self.EL_i) * rand(len(self.I_pop))
        self.setBackgroundInput(self.o.Iext_e*amp, self.o.Iext_i*amp)
        self.E_pop.tau_ad = (self.o.ad_tau_e_mean +
                self.o.ad_tau_e_std*np.random.randn(len(self.E_pop.tau_ad)))*second
        self.E_pop.EL = (self.o.EL_e - self.o.EL_e_spread/2. +
                self.o.EL_e_spread*rand(len(self.E_pop))) * volt
        self.E_pop.taum = (self.o.taum_e - self.o.taum_e_spread/2. +
                self.o.taum_e_spread*rand(len(self.E_pop))) * second
        self.I_pop.tau_ad = (self.o.ad_tau_i_mean +
                self.o.ad_tau_i_std*np.random.randn(len(self.I_pop.tau_ad)))*second
        self.I_pop.EL = (self.o.EL_i - self.o.EL_i_spread/2.0 +
                self.o.EL_i_spread*rand(len(self.I_pop))) * second
        self.I_pop.taum = (self.o.taum_i - self.o.taum_i_spread/2. +
                self.o.taum_i_spread*rand(len(self.I_pop))) * second

    def reinit(self):
        self.net.reinit(states=True)
        self.initStates()

    def _createBasicNet(self, o, clk):

        noise_sigma=o.noise_sigma*volt

        # Setup neuron equations
        # Using exponential integrate and fire model
        # Excitatory population
        Rm_e = o.Rm_e*ohm
        taum_e = o.taum_e*second
        Ce = taum_e/Rm_e
        self.EL_e = o.EL_e*volt
        deltaT_e = o.deltaT_e*volt
        self.Vt_e = o.Vt_e*volt
        Vr_e = o.Vr_e*volt
        g_ahp_e = o.g_ahp_e * siemens
        tau_GABA_rise = o.tau_GABA_rise*second
        tau_GABA_fall = o.tau_GABA_fall*second
        tau1_GABA = tau_GABA_fall
        tau2_GABA = tau_GABA_rise*tau_GABA_fall / (tau_GABA_rise + tau_GABA_fall);
        self.B_GABA = 1/((tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau1_GABA) - 
                (tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau2_GABA))
        
        Vrev_GABA = o.Vrev_GABA*volt


        self.eqs_e = Equations('''
            dvm/dt = 1/C*Im + (noise_sigma*xi/taum_mean**.5): volt
            Im = gL*(EL-vm)*(1+g_ad/gL) + g_ahp*(Eahp - vm) + gL*deltaT*exp((vm-Vt)/deltaT) + Isyn + Iext  : amp
            Isyn = (gi1 - gi2)*(Esyn - vm) : amp
            Iclamp = (gi1 - gi2)*(Esyn - Vclamp) : amp
            dgi1/dt = -gi1/syn_tau1 : siemens
            dgi2/dt = -gi2/syn_tau2 : siemens
            dg_ad/dt = -g_ad/tau_ad : siemens
            dg_ahp/dt = -g_ahp/tau_ahp : siemens
            tau_ad : second
            Iext : amp
            EL : volt
            taum : second
            ''',
            C=Ce,
            gL=1/Rm_e,
            noise_sigma=noise_sigma,
            deltaT=deltaT_e,
            Vt=self.Vt_e,
            Esyn=Vrev_GABA,
            Vclamp = o.Vclamp*volt,
            syn_tau1=tau1_GABA,
            syn_tau2=tau2_GABA,
            taum_mean=o.taum_e*second,
            tau_ahp=o.tau_ahp_e*second,
            Eahp=o.Eahp_e * volt)


        # Inhibitory population
        Rm_i = o.Rm_i*ohm
        taum_i = o.taum_i*second
        Ci = taum_i/Rm_i
        self.EL_i = o.EL_i*volt
        deltaT_i = o.deltaT_i*volt
        self.Vt_i = o.Vt_i*volt
        Vr_i = o.Vr_i*volt
        Vrev_AMPA = o.Vrev_AMPA*volt
        tau_AMPA = o.tau_AMPA*second
        tau_ad_i = o.ad_tau_i_mean*second
        
        self.eqs_i = Equations('''
            dvm/dt = 1/C*Im + (noise_sigma*xi/taum_mean**.5): volt
            Im = gL*(EL-vm)*(1+g_ad/gL) + gL*deltaT*exp((vm-Vt)/deltaT) + Isyn + Iext  : amp
            Isyn = ge*(Esyn - vm) : amp
            Iclamp = ge*(Esyn - Vclamp) : amp
            dge/dt = -ge/syn_tau : siemens
            dg_ad/dt = -g_ad/tau_ad : siemens
            tau_ad : second
            Iext : amp
            EL : volt
            taum : second
            ''',
            C=Ci,
            gL=1/Rm_i,
            noise_sigma=noise_sigma,
            deltaT=deltaT_i,
            Vt=self.Vt_i,
            Esyn=Vrev_AMPA,
            Vclamp=o.Vclamp*volt,
            syn_tau=tau_AMPA,
            taum_mean=o.taum_i*second)


        # Other constants
        refrac_abs = o.refrac_abs*second
        spike_detect_th = o.spike_detect_th*volt



        # Setup neuron groups and connections
        self.E_pop = NeuronGroup(
                N = self.net_Ne,
                model=self.eqs_e,
                threshold=spike_detect_th,
                reset="vm=Vr_e; g_ahp=g_ahp_e",
                refractory=refrac_abs,
                clock=clk)

        self.I_pop = NeuronGroup(
                N = self.net_Ni,
                model=self.eqs_i,
                threshold=spike_detect_th,
                reset=Vr_i,
                refractory=refrac_abs,
                clock=clk)


        self.net = Network(
                self.E_pop,
                self.I_pop)

        # Setup adaptation connections: neuron on itself
        if o.ad_e_g_inc != 0.0:
            self.adaptConn_e = IdentityConnection(self.E_pop, self.E_pop,  'g_ad',
                    weight=o.ad_e_g_inc*siemens)
            self.net.add(self.adaptConn_e)
        if o.ad_i_g_inc != 0.0:
            self.adaptConn_i = IdentityConnection(self.I_pop, self.I_pop, 'g_ad',
                    weight=o.ad_i_g_inc*siemens)
            self.net.add(self.adaptConn_i)



        self.o = o
        self.initStates()


    def setBackgroundInput(self, Iext_e, Iext_i):
        self.E_pop.Iext = linspace(Iext_e, Iext_e, len(self.E_pop))
        self.I_pop.Iext = linspace(Iext_i, Iext_i, len(self.I_pop))

    def _generate_pAMPA_template1D(self, N, mu, sigma):
        '''Generate AMPA probability profile function on an interval [0, N-1]. For
        now it will be a Gaussian profile, with mean mu and std. dev. sigma. The
        boundaries are wrapped-around (N==0). This is distance-dependent
        profile, which take wrap around boundaries into account.'''
        self.pAMPA_templ = np.exp(-(np.arange(N/2+1) - mu)**2/2.0/sigma**2)
        self.pAMPA_templ = np.concatenate((self.pAMPA_templ,
            self.pAMPA_templ[1:N-len(self.pAMPA_templ)+1][::-1]))

    def _generate_pAMPA_template2D(self, N, mu, sigma):
        mg = np.mgrid[0:N/2+1, 0:N/2+1]
        self.pAMPA_templ = np.exp(-(np.sqrt(mg[0]**2 + mg[1]**2) - mu)**2/2.0/sigma**2)

        # Stack horizontally and vertically
        ncols = len(self.pAMPA_templ[0])
        self.pAMPA_templ = np.hstack((self.pAMPA_templ, self.pAMPA_templ[:,
            1:N-ncols+1][:, ::-1]))
        nrows = len(self.pAMPA_templ)
        self.pAMPA_templ = np.vstack((self.pAMPA_templ,
            self.pAMPA_templ[1:N-nrows+1, :][::-1, :]))

    def _generate_pGABA_template1D(self, N, sigma):
        '''Generate GABA probability profile function on an interval [0, N-1].
        It is similar to generate_pAMPA_template but the mean is 0'''
        self.pGABA_templ = np.exp(-(np.arange(N/2+1))**2/2.0/sigma**2)
        self.pGABA_templ = np.concatenate((self.pGABA_templ,
            self.pGABA_templ[1:N-len(self.pGABA_templ)+1][::-1]))


    def _generate_pGABA_template2D(self, N, sigma):
        mg = np.mgrid[0:N/2+1, 0:N/2+1]
        self.pGABA_templ = np.exp(-(mg[0]**2 + mg[1]**2)/2.0/sigma**2)

        # Stack horizontally and vertically
        ncols = len(self.pGABA_templ[0])
        self.pGABA_templ = np.hstack((self.pGABA_templ, self.pGABA_templ[:,
            1:N-ncols+1][:, ::-1]))
        nrows = len(self.pGABA_templ)
        self.pGABA_templ = np.vstack((self.pGABA_templ,
            self.pGABA_templ[1:N-nrows+1, :][::-1, :]))

        
    def connMexicanHat(self, pAMPA_mu, pAMPA_sigma, pGABA_sigma):
        '''Create excitatory and inhibitory connections, Mexican hat ring model.
            pAMPA_mu    Mean of the AMPA Gaussian profile of connection
                        probability
            pAMPA_sigma Std. dev. of the AMPA Gaussian profile of connection
                        probability
            pGABA_sigma Std. dev. of the GABA Gaussian profile of connection
                        probability (Mean of GABA connections is local)
            ndim        Number of dimensions (1 or 2)
        '''

        g_AMPA_mean = self.o.g_AMPA_total/self.net_Ne
        g_AMPA_sigma = np.sqrt(np.log(1 + self.o.g_AMPA_std**2/g_AMPA_mean**2))
        g_AMPA_mu = np.log(g_AMPA_mean) - 1/2*g_AMPA_sigma**2
        g_GABA_mean = self.o.g_GABA_total / self.net_Ni * siemens

        # Generate connection-probability profile functions for GABA and AMPA connections
        self.pAMPA_mu = pAMPA_mu * self.o.Ni
        self.pAMPA_sigma = pAMPA_sigma * self.o.Ni
        self.pGABA_sigma = pGABA_sigma * self.o.Ne

        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge')
        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        if (self.o.ndim ==1):
            self._generate_pAMPA_template1D(self.o.Ni, self.pAMPA_mu, self.pAMPA_sigma)
            self._generate_pGABA_template1D(self.o.Ne, self.pGABA_sigma)
            for i in xrange(self.o.Ne):
                e_norm = int(round(np.double(i)/self.o.Ne*self.o.Ni))

                tmp_templ = np.roll(self.pAMPA_templ, e_norm)
                self.AMPA_conn.W.rows[i] = list((rand(self.o.Ni) <
                    self.o.AMPA_density*tmp_templ).nonzero()[0])
                self.AMPA_conn.W.data[i] = np.random.lognormal(g_AMPA_mu,
                        g_AMPA_sigma, len(self.AMPA_conn.W.rows[i]))*siemens

            for i in xrange(self.o.Ni):
                i_norm = int(round(np.double(i)/self.o.Ni*self.o.Ne))
                

                tmp_templ = np.roll(self.pGABA_templ, i_norm)
                self.GABA_conn1.W.rows[i] = list((rand(self.o.Ne) <
                    self.o.GABA_density*tmp_templ).nonzero()[0])
                self.GABA_conn1.W.data[i] = [self.B_GABA*g_GABA_mean] * len(self.GABA_conn1.W.rows[i])

        elif (self.o.ndim == 2):
            self._generate_pAMPA_template2D(self.o.Ni, self.pAMPA_mu, self.pAMPA_sigma)
            self._generate_pGABA_template2D(self.o.Ne, self.pGABA_sigma)

            # Do the same rolling as in the case of ndim==1, but in two
            # dimensions
            for r in xrange(self.o.Ne):
                r_e_norm = int(round(np.double(r)/self.o.Ne*self.o.Ni))
                tmp_r_templ = np.roll(self.pAMPA_templ, r_e_norm, 0)
                for c in xrange(self.o.Ne):
                    c_e_norm = int(round(np.double(c)/self.o.Ne*self.o.Ni))

                    tmp_templ = np.roll(tmp_r_templ, c_e_norm, 1).ravel('C')
                    it = r*self.o.Ne + c
                    self.AMPA_conn.W.rows[it] = list((rand(self.o.Ni**2) <
                        self.o.AMPA_density*tmp_templ).nonzero()[0])
                    self.AMPA_conn.W.data[it] = np.random.lognormal(g_AMPA_mu,
                            g_AMPA_sigma, len(self.AMPA_conn.W.rows[it]))*siemens

            for r in xrange(self.o.Ni):
                r_i_norm = int(round(np.double(r)/self.o.Ni*self.o.Ne))
                tmp_r_templ = np.roll(self.pGABA_templ, r_i_norm, axis=0)
                for c in xrange(self.o.Ni):
                    c_e_norm = int(round(np.double(c)/self.o.Ni*self.o.Ne))

                    tmp_templ = np.roll(tmp_r_templ, c_e_norm,
                            axis=1).ravel('C')
                    it = r*self.o.Ni + c
                    self.GABA_conn1.W.rows[it] = list((rand(self.o.Ne**2) <
                        self.o.GABA_density*tmp_templ).nonzero()[0])
                    self.GABA_conn1.W.data[it] = [self.B_GABA*g_GABA_mean] * len(self.GABA_conn1.W.rows[it])

        else:
            raise Exception("Number of Mexican hat dimensions must be 0 or 1" +
                ", not" + str(ndim) + ".")


        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)

        self.net.add(self.AMPA_conn, self.GABA_conn1, self.GABA_conn2)
        self.ndim = ndim


    def connRandom(self, AMPA_density, GABA_density):
        '''
        Connect each pair of neurons with a random probability given by:
          
          E --> I given by AMPA_density <0, 1>
          I --> E given by GABA_density <0, 1>
        '''

        g_AMPA_mean = self.o.g_AMPA_total / self.net_Ne / AMPA_density
        g_AMPA_sigma = np.sqrt(np.log(1 + self.o.g_AMPA_std**2/g_AMPA_mean**2))
        g_AMPA_mu = np.log(g_AMPA_mean) - 1/2*g_AMPA_sigma**2
        g_GABA_mean = self.o.g_GABA_total / self.net_Ni / GABA_density * siemens

        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge')
        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        self.AMPA_conn.connect_random(self.E_pop, self.I_pop, AMPA_density, 
                weight=lambda: np.random.lognormal(g_AMPA_mu, g_AMPA_sigma)*siemens)
        self.GABA_conn1.connect_random(self.I_pop, self.E_pop, GABA_density,
                weight=self.B_GABA*g_GABA_mean)
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W) 
  
        self.net.add(self.AMPA_conn, self.GABA_conn1, self.GABA_conn2) 
