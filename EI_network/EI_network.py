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

import logging as lg

lg.basicConfig(level=lg.DEBUG)

# Get a preferred direction for a neuron
def getPreferredDirection(pos_x, pos_y):
# pos_x/y - position of neuron in 2d sheet
    pos4_x = pos_x % 2
    pos2_y = pos_y % 2
    if pos4_x == 0:
        if pos2_y == 0:
            return [-1, 0] # Left
        else:
            return [0, -1] # Down
    else:
        if pos2_y == 0:
            return [0, 1] # up
        else:
            return [1, 0] # Right


class EI_Network:
    def __init__(self, o, clk):
        if o.ndim == 1 or o.ndim == -1:
            self.net_Ne = o.Ne
            self.net_Ni = o.Ni
        elif o.ndim == 2:
            self.net_Ne = o.Ne**2
            self.net_Ni = o.Ni**2
        elif o.ndim == 'twisted_torus':
            self.y_dim = np.sqrt(3)/2.0
            self.Ne_x = o.Ne
            self.Ne_y = int(np.ceil(o.Ne * self.y_dim)) // 2 * 2

            self.Ni_x = o.Ni
            self.Ni_y = int(np.ceil(o.Ni * self.y_dim)) // 2 * 2

            self.net_Ne = self.Ne_x * self.Ne_y
            self.net_Ni = self.Ni_x * self.Ni_y
        else:
            raise Exception("Number of Mexican hat dimensions must be 0 or 1 or"
                    "a 'twisted_torus' string not" + str(ndim) + ".")
        self._createBasicNet(o, clk)

    def initStates(self):
        # Initialize membrane potential randomly
        self.E_pop.vm = self.EL_e + (self.Vt_e-self.EL_e) * rand(len(self.E_pop))
        self.I_pop.vm = self.EL_i + (self.Vt_i-self.EL_i) * rand(len(self.I_pop))
        self.setBackgroundInput(self.o.Iext_e*amp, self.o.Iext_i*amp)
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
            Im = gL*(EL-vm) + g_ahp*(Eahp - vm) + gL*deltaT*exp((vm-Vt)/deltaT) + Isyn + Iext  : amp
            Isyn = (gi1 - gi2)*(Esyn - vm) : amp
            Iclamp = -(gi1 - gi2)*(Esyn - Vclamp) : amp
            dgi1/dt = -gi1/syn_tau1 : siemens
            dgi2/dt = -gi2/syn_tau2 : siemens
            dg_ahp/dt = -g_ahp/tau_ahp : siemens
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
            Isyn = ge*(Esyn - vm) + gNMDA*(Esyn - vm): amp
            Iclamp = -(ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp)): amp
            dge/dt = -ge/syn_tau : siemens
            dg_ad/dt = -g_ad/tau_ad : siemens
            tau_ad : second
            Iext : amp
            EL : volt
            taum : second
            dgNMDA/dt = -gNMDA/(100*msecond) : siemens
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

    def _generate_pAMPA_template2D(self, a, others, pAMPA_mu, pAMPA_sigma, prefDir):
        '''
        !!!
        a ... vector of 2D position of current neuron, the shape must be (1,2)
        b ... matrix of 2D positions of N other neurons, the shape must be (N, 2)
        l_net ... size of the field (see Burak and Fiete, PLoS Comp. B.)
        prefDir ... preferred Direction of a neuron (this will shift the
                    Gaussian profile in the direction of prefDir
        '''
        #gamma = 3.15/l_net**2
        #beta = 3./l_net**2
        ## B is a normalisation constant so max of the profile function is 1
        #B = -((gamma/beta)**(gamma/(beta-gamma)) -
        #        (gamma/beta)**(beta/(beta-gamma)))**-1

        d = a - others #- np.array(prefDir)/self.o.Ni
        d = np.abs(1./2./np.pi *
                np.angle(np.exp(1j*2.*np.pi*d)))
        d2 = d[:, 0]**2 + d[:, 1]**2
        #return -B*(np.exp(-gamma*d2) - np.exp(-beta*d2))
        d = np.sqrt(d2)
        return np.exp(-(d - pAMPA_mu)**2/2/pAMPA_sigma**2)

    def _remap_twisted_torus(self, a, others, prefDir):
        a_x = a[0, 0]
        a_y = a[0, 1]
        prefDir_x = prefDir[0, 0]
        prefDir_y = prefDir[0, 1]

        others_x = others[:, 0] + prefDir_x
        others_y = others[:, 1] + prefDir_y

        d1 = sqrt((a_x - others_x)**2 + (a_y - others_y)**2)
        d2 = sqrt((a_x - others_x - 1.)**2 + (a_y - others_y)**2)
        d3 = sqrt((a_x - others_x + 1.)**2 + (a_y - others_y)**2)
        d4 = sqrt((a_x - others_x + 0.5)**2 + (a_y - others_y - self.y_dim)**2)
        d5 = sqrt((a_x - others_x - 0.5)**2 + (a_y - others_y - self.y_dim)**2)
        d6 = sqrt((a_x - others_x + 0.5)**2 + (a_y - others_y + self.y_dim)**2)
        d7 = sqrt((a_x - others_x - 0.5)**2 + (a_y - others_y + self.y_dim)**2)
        
        return np.min((d1, d2, d3, d4, d5, d6, d7), 0)
            
        # Straight torus here
        #dx = np.min((
        #    np.abs(a_x - others_x),
        #    np.abs(a_x - (others_x + 1.)),
        #    np.abs(a_x - (others_x - 1.))), 0)
        #dy = np.min((
        #    np.abs(a_y - others_y),
        #    np.abs(a_y - (others_y + self.y_dim)),
        #    np.abs(a_y - (others_y - self.y_dim))), 0)
        #return np.sqrt(dx**2 + dy**2)

    def _generate_pAMPA_twisted_torus(self, a, others, pAMPA_mu, pAMPA_sigma,
            prefDir, prefDirC):
        '''
        Here we assume that X coordinates are normalised to <0, 1), and Y
        coordinates are normalised to <0, sqrt(3)/2)
        Y coordinates are twisted, i.e. X will have additional position shifts
        when determining minimum
        '''
        prefDir = np.array(prefDir, dtype=float)
        prefDir[0, 0] = 1. * prefDir[0, 0] * prefDirC / self.Ni_x
        prefDir[0, 1] = 1. * prefDir[0, 1] * prefDirC / self.Ni_y * self.y_dim
        d = self._remap_twisted_torus(a, others, prefDir)
        return np.exp(-(d - pAMPA_mu)**2/2/pAMPA_sigma**2)

    def _generate_pGABA_twisted_torus(self, a, others, sigma, prefDir,
        prefDirC):

        #import pdb; pdb.set_trace()
        prefDir = np.array(prefDir, dtype=float)
        prefDir[0, 0] = 1. * prefDir[0, 0] * prefDirC / self.Ne_x
        prefDir[0, 1] = 1. * prefDir[0, 1] * prefDirC / self.Ne_y * self.y_dim
        d = self._remap_twisted_torus(a, others, prefDir)
        return np.exp(-d**2/2./sigma**2)



    def _generate_pGABA_template2D(self, a, others, sigma, prefDir, prefDirC):
        d = a - others - prefDirC*np.array(prefDir)/self.o.Ne
        d = np.abs(1./2./np.pi *
                np.angle(np.exp(1j*2.*np.pi*d)))
        d2 = d[:, 0]**2 + d[:, 1]**2
        return np.exp(-d2/2./sigma**2)

    def connMexicanHat(self, pAMPA_mu, pAMPA_sigma, pGABA_sigma):
        '''Create excitatory and inhibitory connections, Mexican hat ring model.
            pGABA_sigma Std. dev. of the GABA Gaussian profile of connection
                        probability (Mean of GABA connections is local)
            ndim        Number of dimensions (1 or 2)
        '''

        g_AMPA_mean = self.o.g_AMPA_total/self.net_Ne
        #GABA_density = 9*np.pi*pGABA_sigma**2
        g_GABA_mean = self.o.g_GABA_total / self.net_Ni * siemens

        # Generate connection-probability profile functions for GABA and AMPA connections
        self.pGABA_sigma = pGABA_sigma * self.o.Ne

        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge',
            structure='dense')
        self.NMDA_conn = Connection(self.E_pop, self.I_pop, 'gNMDA',
            structure='dense')
        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        if (self.o.ndim == 2):
            #self.prefDirs = np.ndarray((self.net_Ne, 2))
            X, Y = np.meshgrid(np.arange(self.o.Ni), np.arange(self.o.Ni))
            others_e = 1. * np.vstack((X.ravel(), Y.ravel())).T / self.o.Ni

            E_W = np.asarray(self.AMPA_conn.W)
            for y in xrange(self.o.Ne):
                y_e_norm = float(y) / self.o.Ne

                for x in xrange(self.o.Ne):
                    it = y*self.o.Ne + x

                    x_e_norm = float(x) / self.o.Ne

                    
                    a = np.array([[x_e_norm, y_e_norm]])
                    #pd = getPreferredDirection(x, y)
                    pd = np.array([0, 0])
                    #self.prefDirs[it, :] = pd
                    tmp_templ = self._generate_pAMPA_template2D(a, others_e,
                            pAMPA_mu, pAMPA_sigma, np.array([[pd[0], pd[1]]]))

                    E_W[it, :] = g_AMPA_mean*tmp_templ*siemens

            conn_th = 1e-3
            self.prefDirs_i = np.ndarray((self.net_Ni, 2))
            X, Y = np.meshgrid(np.arange(self.o.Ne), np.arange(self.o.Ne))
            others_i = 1. * np.vstack((X.ravel(), Y.ravel())).T / self.o.Ne
            for y in xrange(self.o.Ni):
                y_i_norm = float(y) / self.o.Ni
                for x in xrange(self.o.Ni):
                    x_i_norm = float(x) / self.o.Ni
                    it = y*self.o.Ni + x

                    a = np.array([[x_i_norm, y_i_norm]])
                    pd = getPreferredDirection(x, y)
                    #pd = np.array([0, 0])
                    self.prefDirs_i[it, :] = pd
                    tmp_templ = self._generate_pGABA_template2D(a, others_i,
                            pGABA_sigma, [[pd[0], pd[1]]], self.o.prefDirC)

                    self.GABA_conn1.W.rows[it] = (tmp_templ >
                            conn_th).nonzero()[0]
                    self.GABA_conn1.W.data[it] = self.B_GABA*g_GABA_mean*tmp_templ[self.GABA_conn1.W.rows[it]]

        elif self.o.ndim == 'twisted_torus':
            X, Y = np.meshgrid(np.arange(self.Ni_x), np.arange(self.Ni_y))
            X = 1. * X / self.Ni_x
            Y = 1. * Y / self.Ni_y * self.y_dim
            others_e = np.vstack((X.ravel(), Y.ravel())).T

            self.prefDirs_e = np.ndarray((self.net_Ne, 2))
            E_W = np.asarray(self.AMPA_conn.W)
            for y in xrange(self.Ne_y):
                y_e_norm = float(y) / self.Ne_y * self.y_dim

                for x in xrange(self.Ne_x):
                    it = y*self.Ne_x + x

                    x_e_norm = float(x) / self.Ne_x

                    
                    a = np.array([[x_e_norm, y_e_norm]])
                    pd = getPreferredDirection(x, y)
                    #pd = np.array([0, 0])
                    self.prefDirs_e[it, :] = pd
                    tmp_templ = self._generate_pAMPA_twisted_torus(a, others_e,
                            pAMPA_mu, pAMPA_sigma, np.array([[pd[0], pd[1]]]),
                            self.o.prefDirC)

                    E_W[it, :] = g_AMPA_mean*tmp_templ*siemens

            conn_th = 1e-5
            self.prefDirs_i = np.ndarray((self.net_Ni, 2))
            X, Y = np.meshgrid(np.arange(self.Ne_x), np.arange(self.Ne_y))
            X = 1. * X / self.Ne_x
            Y = 1. * Y / self.Ne_y * self.y_dim
            others_i = np.vstack((X.ravel(), Y.ravel())).T
            for y in xrange(self.Ni_y):
                y_i_norm = float(y) / self.Ni_y * self.y_dim
                for x in xrange(self.Ni_x):
                    x_i_norm = float(x) / self.Ni_x
                    it = y*self.Ni_x + x

                    a = np.array([[x_i_norm, y_i_norm]])
                    #pd = getPreferredDirection(x, y)
                    pd = np.array([0, 0])
                    self.prefDirs_i[it, :] = pd
                    tmp_templ = self._generate_pGABA_twisted_torus(a, others_i,
                            pGABA_sigma, [[pd[0], pd[1]]], self.o.prefDirC)

                    self.GABA_conn1.W.rows[it] = (tmp_templ >
                            conn_th).nonzero()[0]
                    self.GABA_conn1.W.data[it] = self.B_GABA*g_GABA_mean*tmp_templ[self.GABA_conn1.W.rows[it]]

        else:
            raise Exception("Number of Mexican hat dimensions must be 2" +
                ", not" + str(ndim) + ".")


        self.NMDA_conn.connect(self.E_pop, self.I_pop,
                self.AMPA_conn.W * .01 * self.o.NMDA_amount)
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)

        self.net.add(self.AMPA_conn, self.NMDA_conn, self.GABA_conn1, self.GABA_conn2)
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

    def randomInhibition(self, total_strength, density):
        '''Random inhibitory connections from I to E only'''

        g_GABA_mean = total_strength / self.net_Ni / density * siemens

        self.extraGABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.extraGABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        self.extraGABA_conn1.connect_random(self.I_pop, self.E_pop, density,
                weight=self.B_GABA*g_GABA_mean)
        self.extraGABA_conn2.connect(self.I_pop, self.E_pop, self.extraGABA_conn1.W) 

        self.net.add(self.extraGABA_conn1, self.extraGABA_conn2)
