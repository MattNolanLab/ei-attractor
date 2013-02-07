#
#   simulation_basic_grids.py
#
#   Main simulation run: grid fields with theta input and all the inhibition
#   (for gamma) and place input.
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

from matplotlib.backends.backend_pdf import PdfPages

from scipy      import linspace
from scipy.io   import loadmat
from scipy.io   import savemat
from optparse   import OptionParser

from parameters              import *
from grid_cell_network_brian import *
from custombrian             import *

import time
import math
import sys
import numpy as np
import logging as lg

import default_params


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--ndumps",  type="int",    help="Number of data output dumps during the simulation")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)



################################################################################
#            Helper functions to generate ring-like and Gaussian
#                  synaptic weights on a twisted torus
################################################################################
def _remap_twisted_torus(self, a, others, prefDir):
    '''
    Take a, the position of one neuron on the twisted torus and compute the
    distance of this neuron from others, accounting from preferred direction of
    the neuron (a). Return an array of distances.
    '''
    a_x = a[0, 0]
    a_y = a[0, 1]
    prefDir_x = prefDir[0, 0]
    prefDir_y = prefDir[0, 1]

    others_x = others[:, 0] + prefDir_x
    others_y = others[:, 1] + prefDir_y

    d1 = sqrt((a_x - others_x)**2 + (a_y - others_y)**2)
    d2 = sqrt((a_x - others_x - 1.)**2 + (a_y - others_y)**2)
    d3 = sqrt((a_x - others_x + 1.)**2 + (a_y - others_y)**2)
    d4 = sqrt((a_x - others_x + 0.5)**2 + (a_y - others_y - y_dim)**2)
    d5 = sqrt((a_x - others_x - 0.5)**2 + (a_y - others_y - y_dim)**2)
    d6 = sqrt((a_x - others_x + 0.5)**2 + (a_y - others_y + y_dim)**2)
    d7 = sqrt((a_x - others_x - 0.5)**2 + (a_y - others_y + y_dim)**2)
    
    return np.min((d1, d2, d3, d4, d5, d6, d7), 0)
        

def _generateRinglikeWeights(self, a, others, mu, sigma, prefDir, prefDirC):
    '''
    Return connection weights given by a ring-like function.

    Here we assume that X coordinates are normalised to <0, 1), and Y
    coordinates are normalised to <0, sqrt(3)/2)
    Y coordinates are twisted, i.e. X will have additional position shifts
    when determining minimum.
    '''
    prefDir = np.array(prefDir, dtype=float) * prefDirC
    d = _remap_twisted_torus(a, others, prefDir)
    return np.exp(-(d - mu)**2/2/sigma**2)


def _generateGaussianWeights(self, a, others, sigma, prefDir, prefDirC):
    '''
    Return connection weights given by a ring-like function.

    Here we assume that X coordinates are normalised to <0, 1), and Y
    coordinates are normalised to <0, sqrt(3)/2)
    Y coordinates are twisted, i.e. X will have additional position shifts
    when determining minimum.
    '''

    #import pdb; pdb.set_trace()
    prefDir = np.array(prefDir, dtype=float) * prefDirC
    d = _remap_twisted_torus(a, others, prefDir)
    return np.exp(-d**2/2./sigma**2)



def _centerSurroundConnection(self, AMPA_gaussian, pAMPA_mu, pAMPA_sigma, pGABA_mu, pGABA_sigma):
    '''
    Create a center-surround excitatory and inhibitory connections between
    both populations.

    The connections are remapped to [1.0, sqrt(3)/2], whether the topology
    is a twisted torus or just a regular torus.

    AMPA_gaussian switches between two cases:
        true    Each exciatory neuron has a 2D excitatory gaussian profile,
                while each inhibitory neuron has a ring-like profile
                pAMPA_mu, pAMPA_sigma, pGABA_sigma are used,
                pGABA_mu is discarded
        false   Each excitatory neuron has a ring-like profile, while
                each inhibitory neuron has a gaussian profile.
                pAMPA_sigma, pGABA_mu, pGABA_sigma are used,
                pAMPA_mu is discarded
    '''
    g_AMPA_mean = no.g_AMPA_total / net_Ne
    g_GABA_mean = no.g_GABA_total / net_Ni

    # E --> I connections
    X, Y = np.meshgrid(np.arange(Ni_x), np.arange(Ni_y))
    X = 1. * X / Ni_x
    Y = 1. * Y / Ni_y * y_dim
    others_e = np.vstack((X.ravel(), Y.ravel())).T

    prefDirs_e = np.ndarray((net_Ne, 2))
    for y in xrange(Ne_y):
        y_e_norm = float(y) / Ne_y * y_dim

        for x in xrange(Ne_x):
            it = y*Ne_x + x

            x_e_norm = float(x) / Ne_x

            a = np.array([[x_e_norm, y_e_norm]])
            pd_e = getPreferredDirection(x, y)
            prefDirs_e[it, :] = pd_e

            pd_norm_e = np.ndarray((1, 2)) # normalise preferred dirs before sending them
            pd_norm_e[0, 0] = 1. * pd_e[0] / Ni_x
            pd_norm_e[0, 1] = 1. * pd_e[1] / Ni_y * y_dim
            if AMPA_gaussian == 1:
                tmp_templ = _generateGaussianWeights(a, others_e,
                        pAMPA_sigma, pd_norm_e, no.prefDirC_e)
            elif AMPA_gaussian == 0:
                tmp_templ = _generateRinglikeWeights(a, others_e,
                        pAMPA_mu, pAMPA_sigma, pd_norm_e, no.prefDirC_e)
            else:
                raise Exception('AMPA_gaussian parameters must be 0 or 1')

            tmp_templ *= g_AMPA_mean
            # tmp_templ down here must be in the proper units (e.g. nS)
            tmp_templ = _addToConnections(tmp_templ, no.condAddPercSynapses_e, no.condAdd_e)
            _divergentConnectEI(it, range(net_Ni), tmp_templ)

    # I --> E connections
    conn_th = 1e-5
    prefDirs_i = np.ndarray((net_Ni, 2))
    X, Y = np.meshgrid(np.arange(Ne_x), np.arange(Ne_y))
    X = 1. * X / Ne_x
    Y = 1. * Y / Ne_y * y_dim
    others_i = np.vstack((X.ravel(), Y.ravel())).T
    for y in xrange(Ni_y):
        y_i_norm = float(y) / Ni_y * y_dim
        for x in xrange(Ni_x):
            x_i_norm = float(x) / Ni_x
            it = y*Ni_x + x

            a = np.array([[x_i_norm, y_i_norm]])
            pd_i = getPreferredDirection(x, y)
            prefDirs_i[it, :] = pd_i

            pd_norm_i = np.ndarray((1, 2)) # normalise preferred dirs before sending them
            pd_norm_i[0, 0] = 1. * pd_i[0] / Ne_x
            pd_norm_i[0, 1] = 1. * pd_i[1] / Ne_y * y_dim
            if AMPA_gaussian == 1:
                tmp_templ = _generateRinglikeWeights(a, others_i,
                        pGABA_mu, pGABA_sigma, pd_norm_i, no.prefDirC_i)
            elif AMPA_gaussian == 0:
                tmp_templ = _generateGaussianWeights(a, others_i,
                        pGABA_sigma, pd_norm_i, no.prefDirC_i)
            else:
                raise Exception('AMPA_gaussian parameters must be 0 or 1')

            E_nid = (tmp_templ > conn_th).nonzero()[0]
            _divergentConnectIE(it, E_nid, B_GABA*g_GABA_mean*tmp_templ[E_nid])




################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)




# Create network
tau1_GABA = no.tau_GABA_A_fall
tau2_GABA = no.tau_GABA_A_rise * no.tau_GABA_A_fall / \
        (no.tau_GABA_A_rise + no.tau_GABA_A_fall);
B_GABA = 1/((tau2_GABA/tau1_GABA)**(no.tau_GABA_A_rise/tau1_GABA) - 
        (tau2_GABA/tau1_GABA)**(no.tau_GABA_A_rise/tau2_GABA))

# Stellate cell equations
eqs_e = Equations('''
    dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)                      : volt
    Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                                : amp
    Im          = gL*(EL-vm) + g_ahp*(Eahp - vm) + Ispike + Isyn + Iext        : amp
    Isyn        = B_GABA*(gi1 - gi2)*(Esyn_i - vm) + ge*(Esyn_e - vm) + gNMDA*(Esyn_e - vm) : amp
    Iclamp      = -(B_GABA*(gi1 - gi2)*(Esyn_i - Vclamp) + ge*(Esyn_e - Vclamp) + gNMDA*(Esyn_e - Vclamp))                        : amp
    Iclamp_all  = -( gL*(EL-Vclamp) + gL*deltaT*exp((Vclamp-Vt)/deltaT) + Iext ) + Iclamp : amp
    dge/dt      = -ge/syn_tau_e                                                : siemens
    dgNMDA/dt   = -gNMDA/tau_NMDA_fall                                         : siemens
    dgi1/dt     = -gi1/syn_tau1                                                : siemens
    dgi2/dt     = -gi2/syn_tau2                                                : siemens
    dg_ahp/dt   = -g_ahp/tau_ahp                                               : siemens
    Iext        = Iext_const + Iext_theta + Iext_vel + Iext_start + Iext_place : amp
    Iext_const                                                                 : amp
    Iext_theta                                                                 : amp
    Iext_vel                                                                   : amp
    Iext_start                                                                 : amp
    Iext_place                                                                 : amp
    EL                                                                         : volt
    taum                                                                       : second
    ''',
    C             = no.taum_e * no.gL_e * pF,
    gL            = no.gL_e * nS,
    noise_sigma   = no.noise_sigma * mV,
    deltaT        = no.deltaT_e * mV,
    Vt            = no.Vt_e * mV,
    Esyn_i        = no.E_GABA_A * mV,
    Esyn_e        = no.E_AMPA * mV,
    Vclamp        = no.Vclamp * mV,
    syn_tau_e     = no.tau_AMPA * ms,
    tau_NMDA_fall = no.tau_NMDA_fall * ms,
    syn_tau1      = tau1_GABA * ms,
    syn_tau2      = tau2_GABA * ms,
    B_GABA        = B_GABA,
    taum_mean     = no.taum_e * ms,
    tau_ahp       = no.tau_AHP_e * ms,
    Eahp          = no.E_AHP_e * mV)


# Interneuron equations
eqs_i = Equations('''
    dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)           : volt
    Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                     : amp
    Im          = gL*(EL-vm)*(1+g_ad/gL) + Ispike + Isyn + Iext     : amp
    Isyn        = ge*(Esyn - vm) + gNMDA*(Esyn - vm)                : amp
    Iclamp      = -(ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp))     : amp
    Iclamp_all  = -( gL*(EL-Vclamp) + gL*deltaT*exp((Vclamp-Vt)/deltaT) + ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp) + Iext )     : amp
    dge/dt      = -ge/syn_tau                                       : siemens
    dg_ad/dt    = -g_ad/tau_ad                                      : siemens
    dgNMDA/dt   = -gNMDA/tau_NMDA_fall                              : siemens
    tau_ad                                                          : second
    Iext        = Iext_const + Iext_theta + Iext_vel                : amp
    Iext_const                                                      : amp
    Iext_theta                                                      : amp
    Iext_vel                                                        : amp
    Iext_start                                                      : amp
    EL                                                              : volt
    taum                                                            : second
    ''',
    C             = no.taum_i * no.gL_i * pF,
    gL            = no.gL_i * nS,
    noise_sigma   = no.noise_sigma * mV,
    deltaT        = no.deltaT_i * mV,
    Vt            = no.Vt_i * mV,
    Esyn          = no.E_AMPA * mV,
    Vclamp        = no.Vclamp * mV,
    syn_tau       = no.tau_AMPA * ms,
    tau_NMDA_fall = no.tau_NMDA_fall * ms,
    taum_mean     = no.taum_i * ms)


# Other constants
g_AHP_e = no.g_AHP_e_max * nS
Vr_e    = no.Vr_e * mV    


# Setup neuron groups and connections
E_pop = NeuronGroup(
        N = net_Ne,
        model=eqs_e,
        threshold=no.V_peak_e * mV,
        reset="vm=Vr_e; g_ahp=g_AHP_e",
        refractory=no.t_ref_e * msecond,
        clock=_simulationClock)

I_pop = NeuronGroup(
        N = net_Ni,
        model=eqs_i,
        threshold=no.V_peak_i * mV,
        reset=no.Vr_i * mV,
        refractory=no.t_ref_i * msecond,
        clock=_simulationClock)

net = Network(E_pop, I_pop)

# Setup adaptation connections: neuron on itself
if no.ad_i_g_inc != 0.0:
    adaptConn_i = IdentityConnection(I_pop, I_pop, 'g_ad',
            weight=no.ad_i_g_inc*nS)
    net.add(adaptConn_i)

# Connect E-->I and I-->E
AMPA_conn = Connection(E_pop, I_pop, 'ge',
    structure='dense')
NMDA_conn = Connection(E_pop, I_pop, 'gNMDA',
    structure='dense')
GABA_conn1 = Connection(I_pop, E_pop, 'gi1')
GABA_conn2 = Connection(I_pop, E_pop, 'gi2')

# Weight matrices which are used in _divergentConnectXY() functions
_E_W = np.asarray(AMPA_conn.W)

_centerSurroundConnection(no.AMPA_gaussian, no.pAMPA_mu,
        no.pAMPA_sigma, no.pGABA_mu, no.pGABA_sigma)

# Now simply copy AMPA --> NMDA and GABA_conn1 --> GABA_conn2
NMDA_conn.connect(E_pop, I_pop, AMPA_conn.W * .01 * no.NMDA_amount)
GABA_conn2.connect(I_pop, E_pop, GABA_conn1.W)

net.add(AMPA_conn, NMDA_conn, GABA_conn1, GABA_conn2)


# Initialisation of neuron states
# Initialize membrane potential randomly
E_pop.vm       = (no.EL_e + (no.Vt_e-no.EL_e) * np.random.rand(len(E_pop))) * mV
E_pop.gi1      = 0.0
E_pop.gi2      = 0.0
E_pop.g_ahp    = 0.0

I_pop.vm       = (no.EL_i + (no.Vt_i-no.EL_i) * np.random.rand(len(I_pop))) * mV
I_pop.ge       = 0.0
I_pop.g_ad     = 0.0
I_pop.gNMDA    = 0.0

# Initialise cellular properties
E_pop.EL       = uniformDistrib(no.EL_e,   no.EL_e_spread,   len(E_pop)) * mV
E_pop.taum     = uniformDistrib(no.taum_e, no.taum_e_spread, len(E_pop)) * ms
I_pop.EL       = uniformDistrib(no.EL_i,   no.EL_i_spread,   len(I_pop)) * mV
I_pop.taum     = uniformDistrib(no.taum_i, no.taum_i_spread, len(I_pop)) * ms

I_pop.tau_ad   = (no.ad_tau_i_mean + no.ad_tau_i_std * np.random.randn(len(I_pop.tau_ad))) * ms


exit(1)


ei_net.setConstantCurrent()
ei_net.setStartCurrent()

ei_net.uniformInhibition()
ei_net.setThetaCurrentStimulation()



#const_v = [0.0, 1.0]
#ei_net.setConstantVelocityCurrent_e(const_v)
ei_net.setVelocityCurrentInput_e()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop)

stateMon_e          = RecentStateMonitor(ei_net.E_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_i          = RecentStateMonitor(ei_net.I_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e   = RecentStateMonitor(ei_net.E_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = RecentStateMonitor(ei_net.I_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = RecentStateMonitor(ei_net.E_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = RecentStateMonitor(ei_net.I_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)


#x_lim = [options.time-0.5, options.time]
x_lim = [options.time/1e3 - 1, options.time/1e3]

################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()
ei_net.net.run(options.time/options.ndumps*msecond, report='stdout')
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


output_fname = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}".format(options.output_dir,
        options.fileNamePrefix, options.job_num, trial_it, dump_it)


F_tstart = 0
F_tend = options.time*1e-3
F_dt = 0.02
F_winLen = 0.25
Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)


# plot firing rates
figure(figsize=figSize)
subplot(211)
T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
pcolormesh(T, FR, Fe)
ylabel('E Neuron no.')
colorbar()
subplot(212)
T, FR = np.meshgrid(Fi_t, np.arange(ei_net.net_Ni))
pcolormesh(T, FR, Fi)
xlabel('Time (s)')
ylabel('I Neuron no.')
colorbar()
savefig(output_fname + '_firing_rate_e.png')

figure()
ax = subplot(211)
plot(stateMon_e.times, stateMon_e.values[:, 0:2]/mV)
ylabel('E membrane potential (mV)')
subplot(212, sharex=ax)
plot(stateMon_i.times, stateMon_i.values[:, 0:2]/mV)
xlabel('Time (s)')
ylabel('I membrane potential (mV)')
xlim(x_lim)
savefig(output_fname + '_Vm.pdf')


figure()
ax = subplot(211)
plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[:, 0:2]/pA)
ylabel('E synaptic current (pA)')
ylim([0, 3000])
subplot(212, sharex=ax)
plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[:, 0:2]/pA)
ylim([-1500, 0])
xlabel('Time (s)')
ylabel('I synaptic current (pA)')
xlim(x_lim)
savefig(output_fname + '_Isyn.pdf')

figure()
ax = subplot(211)
plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[:, 1]/pA)
ylabel('E external current (pA)')
subplot(212, sharex=ax)
plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[:, 0]/pA)
xlabel('Time (s)')
ylabel('I external current (pA)')
xlim(x_lim)
savefig(output_fname + '_Iext.pdf')

figure()
pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
xlabel('E neuron no.')
ylabel('E neuron no.')
colorbar()
axis('equal')
savefig(output_fname + '_firing_snapshot_e.png')

figure()
pcolormesh(np.reshape(Fi[:, len(Fi_t)/2], (ei_net.Ni_y, ei_net.Ni_x)))
xlabel('I neuron no.')
ylabel('I neuron no.')
colorbar()
axis('equal')
savefig(output_fname + '_firing_snapshot_i.png')

# Print a plot of bump position
F_dt = 0.02
F_winLen = 0.25
(pos, bumpPos_times) = spikeMon_e.torusPopulationVector([ei_net.Ne_x,
    ei_net.Ne_y], options.theta_start_t*1e-3, options.time*1e-3, F_dt, F_winLen)
figure(figsize=figSize)
plot(bumpPos_times, pos)
xlabel('Time (s)')
ylabel('Bump position (neurons)')
ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])

savefig(output_fname + '_bump_position.pdf')

        

#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

