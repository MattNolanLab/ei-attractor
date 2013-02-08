#
#   gridcells_brian_demo.py
#
#   Grid cell network demonstration script.
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
from brian      import *

import time
import numpy as np

import parameters


################################################################################
#                      Options and global settings
################################################################################
opt = parameters.defaultParameters

# Setup neuron numbers for each dimension (X, Y)
# We have a single bump and to get hexagonal receptive fields the X:Y
# size ratio must be 1:sqrt(3)/2
y_dim = np.sqrt(3)/2.0
Ne_x = opt.Ne
Ne_y = int(np.ceil(opt.Ne * y_dim)) // 2 * 2

Ni_x = opt.Ni
Ni_y = int(np.ceil(opt.Ni * y_dim)) // 2 * 2

net_Ne = Ne_x * Ne_y
net_Ni = Ni_x * Ni_y


# Clocks
simulationClock = Clock(dt = opt.sim_dt * msecond)


# Other
figSize = (12,8)



################################################################################
#            Helper functions to generate ring-like and Gaussian
#                  synaptic weights on a twisted torus
################################################################################
def uniformDistrib(mean, spread, N):
    return mean - spread/2.0 * rand(N)


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


def _remap_twisted_torus(a, others, prefDir):
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
        

def _generateRinglikeWeights(a, others, mu, sigma, prefDir, prefDirC):
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


def _generateGaussianWeights(a, others, sigma, prefDir, prefDirC):
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



def _centerSurroundConnection(pAMPA_mu, pAMPA_sigma, pGABA_sigma, AMPA_conn, GABA_conn):
    '''
    Create a center-surround excitatory and inhibitory connections between
    both populations.

    The connections are remapped to [1.0, sqrt(3)/2], whether the topology
    is a twisted torus or just a regular torus.
    '''
    g_AMPA_mean = float(opt.g_AMPA_total) / net_Ne
    g_GABA_mean = float(opt.g_GABA_total) / net_Ni

    # E --> I connections
    X, Y = np.meshgrid(np.arange(Ni_x), np.arange(Ni_y))
    X = 1. * X / Ni_x
    Y = 1. * Y / Ni_y * y_dim
    others_e = np.vstack((X.ravel(), Y.ravel())).T

    _E_W = np.asarray(AMPA_conn.W) # Direct setting of connections into a dense array

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
            tmp_templ = _generateRinglikeWeights(a, others_e,
                    pAMPA_mu, pAMPA_sigma, pd_norm_e, opt.prefDirC_e)

            tmp_templ *= g_AMPA_mean
            # tmp_templ down here must be in the proper units (e.g. nS)
            # _E_W[X, Y] --> X...pre, Y...postsynaptic
            _E_W[it, range(net_Ni)] = tmp_templ * nS


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
            tmp_templ = _generateGaussianWeights(a, others_i,
                    pGABA_sigma, pd_norm_i, opt.prefDirC_i)

            E_nid = (tmp_templ > conn_th).nonzero()[0]
            GABA_conn.W.rows[it] = E_nid
            GABA_conn.W.data[it] = B_GABA*g_GABA_mean*tmp_templ[E_nid] * nS


def spikeMonitorFiringRate(spikeMonitor, tstart, tend):
    N = len(spikeMonitor.source)
    firingRate = np.ndarray((N, ))
    for n_it in xrange(N):
        spikes = spikeMonitor[n_it]
        firingRate[n_it] = np.sum(np.logical_and(spikes >= tstart, spikes <= tend)) / (tend - tstart)

    return firingRate



################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()


# Create network
tau1_GABA = opt.tau_GABA_A_fall
tau2_GABA = opt.tau_GABA_A_rise * opt.tau_GABA_A_fall / \
        (opt.tau_GABA_A_rise + opt.tau_GABA_A_fall);
B_GABA = 1/((tau2_GABA/tau1_GABA)**(opt.tau_GABA_A_rise/tau1_GABA) - 
        (tau2_GABA/tau1_GABA)**(opt.tau_GABA_A_rise/tau2_GABA))

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
    C             = opt.taum_e * opt.gL_e * pF,
    gL            = opt.gL_e * nS,
    noise_sigma   = opt.noise_sigma * mV,
    deltaT        = opt.deltaT_e * mV,
    Vt            = opt.Vt_e * mV,
    Esyn_i        = opt.E_GABA_A * mV,
    Esyn_e        = opt.E_AMPA * mV,
    Vclamp        = opt.Vclamp * mV,
    syn_tau_e     = opt.tau_AMPA * ms,
    tau_NMDA_fall = opt.tau_NMDA_fall * ms,
    syn_tau1      = tau1_GABA * ms,
    syn_tau2      = tau2_GABA * ms,
    B_GABA        = B_GABA,
    taum_mean     = opt.taum_e * ms,
    tau_ahp       = opt.tau_AHP_e * ms,
    Eahp          = opt.E_AHP_e * mV)


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
    C             = opt.taum_i * opt.gL_i * pF,
    gL            = opt.gL_i * nS,
    noise_sigma   = opt.noise_sigma * mV,
    deltaT        = opt.deltaT_i * mV,
    Vt            = opt.Vt_i * mV,
    Esyn          = opt.E_AMPA * mV,
    Vclamp        = opt.Vclamp * mV,
    syn_tau       = opt.tau_AMPA * ms,
    tau_NMDA_fall = opt.tau_NMDA_fall * ms,
    taum_mean     = opt.taum_i * ms)


# Other constants
g_AHP_e = opt.g_AHP_e_max * nS
Vr_e    = opt.Vr_e * mV    


# Setup neuron groups and connections
E_pop = NeuronGroup(
        N = net_Ne,
        model=eqs_e,
        threshold=opt.V_peak_e * mV,
        reset="vm=Vr_e; g_ahp=g_AHP_e",
        refractory=opt.t_ref_e * msecond,
        clock=simulationClock)

I_pop = NeuronGroup(
        N = net_Ni,
        model=eqs_i,
        threshold=opt.V_peak_i * mV,
        reset=opt.Vr_i * mV,
        refractory=opt.t_ref_i * msecond,
        clock=simulationClock)

net = Network(E_pop, I_pop)

# Setup adaptation connections: neuron on itself
if opt.ad_i_g_inc != 0.0:
    adaptConn_i = IdentityConnection(I_pop, I_pop, 'g_ad', weight=opt.ad_i_g_inc*nS)
    net.add(adaptConn_i)

# Connect E-->I and I-->E
AMPA_conn = Connection(E_pop, I_pop, 'ge', structure='dense')
NMDA_conn = Connection(E_pop, I_pop, 'gNMDA', structure='dense')
GABA_conn1 = Connection(I_pop, E_pop, 'gi1')
GABA_conn2 = Connection(I_pop, E_pop, 'gi2')

_centerSurroundConnection(opt.pAMPA_mu, opt.pAMPA_sigma, opt.pGABA_sigma, AMPA_conn, GABA_conn1)

# Now simply copy AMPA --> NMDA and GABA_conn1 --> GABA_conn2
NMDA_conn.connect(E_pop, I_pop, AMPA_conn.W * .01 * opt.NMDA_amount)
GABA_conn2.connect(I_pop, E_pop, GABA_conn1.W)

net.add(AMPA_conn, NMDA_conn, GABA_conn1, GABA_conn2)


# Initialisation of neuron states
# Initialize membrane potential randomly
E_pop.vm       = (opt.EL_e + (opt.Vt_e-opt.EL_e) * np.random.rand(len(E_pop))) * mV
E_pop.gi1      = 0.0
E_pop.gi2      = 0.0
E_pop.g_ahp    = 0.0

I_pop.vm       = (opt.EL_i + (opt.Vt_i-opt.EL_i) * np.random.rand(len(I_pop))) * mV
I_pop.ge       = 0.0
I_pop.g_ad     = 0.0
I_pop.gNMDA    = 0.0

# Initialise cellular properties
E_pop.EL       = uniformDistrib(opt.EL_e,   opt.EL_e_spread,   len(E_pop)) * mV
E_pop.taum     = uniformDistrib(opt.taum_e, opt.taum_e_spread, len(E_pop)) * ms
I_pop.EL       = uniformDistrib(opt.EL_i,   opt.EL_i_spread,   len(I_pop)) * mV
I_pop.taum     = uniformDistrib(opt.taum_i, opt.taum_i_spread, len(I_pop)) * ms

I_pop.tau_ad   = (opt.ad_tau_i_mean + opt.ad_tau_i_std * np.random.randn(len(I_pop.tau_ad))) * ms


################################################################################
#                      External input currents settings
################################################################################
E_pop.Iext_const = opt.Iext_e_const * pA
I_pop.Iext_const = opt.Iext_i_const * pA


# Bump initialisation currents
_startCurrentClock = Clock(dt=50*ms)

init_X, init_Y = np.meshgrid(np.arange(Ne_x), np.arange(Ne_y))
init_X = init_X - Ne_x/2.0
init_Y = init_Y - Ne_y/2.0
Istart_coeff = np.array(np.sqrt(init_X**2 + init_Y**2) < opt.Iext_start_size*Ne_y, dtype=float).flatten()

@network_operation(_startCurrentClock)
def startCurrentFun():
    if simulationClock.t < opt.Iext_start_dur*msecond:
        E_pop.Iext_start = Istart_coeff * opt.Iext_start * pA
        print "Bump initialisation..."
    else:
        E_pop.Iext_start = 0.0

net.add(startCurrentFun)



#const_v = [0.0, 1.0]
#setConstantVelocityCurrent_e(const_v)
#setVelocityCurrentInput_e()


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################


state_record_e = [Ne_x/2 -1 , Ne_y/2*Ne_x + Ne_x/2 - 1]
state_record_i = [Ni_x/2 - 1, Ni_y/2*Ni_x + Ni_x/2 - 1]

spikeMon_e          = SpikeMonitor(E_pop)
spikeMon_i          = SpikeMonitor(I_pop)

stateMon_e          = RecentStateMonitor(E_pop, 'vm',     duration=opt.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_i          = RecentStateMonitor(I_pop, 'vm',     duration=opt.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e   = RecentStateMonitor(E_pop, 'Iclamp', duration=opt.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = RecentStateMonitor(I_pop, 'Iclamp', duration=opt.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = RecentStateMonitor(E_pop, 'Iext',   duration=opt.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = RecentStateMonitor(I_pop, 'Iext',   duration=opt.stateMonDur*ms,   record = state_record_i, clock=simulationClock)

net.add(spikeMon_e, spikeMon_i)
net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
net.add(stateMon_Iext_e, stateMon_Iext_i)


#x_lim = [opt.time-0.5, opt.time]
x_lim = [opt.time/1e3 - 1, opt.time/1e3]

################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()
net.run(opt.time*msecond, report='stdout')
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


output_fname = "{0}/gridcells_sim".format(opt.output_dir)



# plot firing rates
figure(figsize=figSize)
ax = subplot(211)
raster_plot(spikeMon_e)
ylabel('E Neuron no.')
subplot(212, sharex=ax)
raster_plot(spikeMon_i)
xlabel('Time (s)')
ylabel('I Neuron no.')
savefig(output_fname + '_raster_plot.png')

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


#figure()
#pcolormesh(np.reshape(AMPA_conn.W[0, :], (Ni_y, Ni_x)))
#colorbar()
#xlabel('I neuron #')
#ylabel('I neuron #')


figure()
ax = subplot(211)
plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[:, 0:2]/pA)
ylabel('E synaptic current (pA)')
ylim([0, 1500])
subplot(212, sharex=ax)
plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[:, 0:2]/pA)
ylim([-1500, 0])
xlabel('Time (s)')
ylabel('I synaptic current (pA)')
xlim(x_lim)
savefig(output_fname + '_Isyn.pdf')


F_tstart = opt.time*1e-3 - 1
F_tend = opt.time*1e-3


figure()
Fe = spikeMonitorFiringRate(spikeMon_e, F_tstart, F_tend)
pcolormesh(np.reshape(Fe, (Ne_y, Ne_x)))
xlabel('E neuron no.')
ylabel('E neuron no.')
colorbar()
axis('equal')
savefig(output_fname + '_firing_snapshot_e.png')

figure()
Fi = spikeMonitorFiringRate(spikeMon_i, F_tstart, F_tend)
pcolormesh(np.reshape(Fi, (Ni_y, Ni_x)))
xlabel('I neuron no.')
ylabel('I neuron no.')
colorbar()
axis('equal')
savefig(output_fname + '_firing_snapshot_i.png')

        

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

