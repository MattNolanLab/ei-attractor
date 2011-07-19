from brian import *
import numpy as np

# Number of neurons
Ne = 800
Ni = 200
N = Ne + Ni

# membrane time constant
taum_e = 10 * ms              
taum_i = 2 * ms

taue = 1 * ms               # excitatory synaptic time constant
taui = 5 * ms              # inhibitory synaptic time constant
Vt = -59 * mV               # spike threshold
Vr = -62 * mV               # reset value
El = -60 * mV               # resting potential
we = (50.0 / N) * mV  # excitatory synaptic weight
wi = (25.0 / N) * mV   # inhibitory synaptic weight

e_sparseness = 0.75
i_sparseness = 0.75

I_e = 1.1 * mV
I_i = 0 * mV

noise_sigma = 0.0 * mV

E_eqs = Equations('''
    dV/dt  = (-gi-(V-El) + I_e)/taum_e + noise_sigma*xi/taum_e**.5 : volt
    dgi/dt = -gi/taui            : volt
    ''')

I_eqs = Equations('''
    dV/dt  = (ge-(V-El) + I_i)/taum_i + noise_sigma*xi/taum_i**.5 : volt
    dge/dt = -ge/taue            : volt
    ''')

#srand(7)

Ge = NeuronGroup(Ne, model=E_eqs, threshold=Vt, reset=Vr)
Gi = NeuronGroup(Ni, model=I_eqs, threshold=Vt, reset=Vr)


Ce = Connection(Ge, Gi, 'ge', sparseness=e_sparseness, weight=we)
Ci = Connection(Gi, Ge, 'gi', sparseness=i_sparseness, weight=wi)

Me = SpikeMonitor(Ge)
Mi = SpikeMonitor(Gi)

MVe = StateMonitor(Ge, 'V', record=0)
MVi = StateMonitor(Gi, 'V', record=0)

Ge.V = Vr + (Vt - Vr) * rand(len(Ge))
Gi.V = Vr + (Vt - Vr) * rand(len(Gi))


###############################################################################
print('Running simulation...')
run(1000*ms)

subplot(411)
raster_plot(Me, title='Pyramidal population', newfigure=False)
subplot(412)
raster_plot(Mi, title='Interneuron population', newfigure=False)

subplot(413)
plot(MVe.times/ms, MVe[0]/mV)
title('Pyramidal neuron')
xlabel('Time (ms)')
ylabel('V (mV)')

subplot(414)
plot(MVi.times/ms, MVi[0]/mV)
xlabel('Time (ms)')
ylabel('V (mV)')


show()
