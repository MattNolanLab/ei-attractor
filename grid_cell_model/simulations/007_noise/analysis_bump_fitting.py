#
#   analysis_bump_fitting.py
#
#   Analysis of fitting the bump and theta/gamma frequency.
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
import numpy as np
from matplotlib.pyplot  import *

from matplotlib.mlab    import detrend_mean, psd, window_hanning

from analysis.spikes    import slidingFiringRateTuple, torusPopulationVector
from analysis.image     import Position2D, fitGaussianBumpTT
from analysis.signal    import relativePower, maxPowerFrequency
from data_storage       import DataStorage



## Compute power of gamma relative to the total power in the signal
#
# @param stateMon   NEST state monitor(s).
# @param gRange     A tuple containing the gamma range (Hz)
# @param sigName    A name of the signal to extract from the state monitor.
# @param tstart     When to start analysis (s)
#
def relativeGammaPower(stateMon, gRange, sigName, tstart=0.0):

    N = len(stateMon)
    stat = stateMon

    relP = np.ndarray((N, ))
    maxF = np.ndarray((N, ))
    powerSpectra_P = []
    powerSpectra_F = []

    for nidx in xrange(N):
        print "relativeGammaPower: ", nidx
        times = stat[nidx]['events']['times'] 
        dt = (times[1] - times[0]) * 1e-3 # sec
        sig = stat[nidx]['events'][sigName][int(tstart/dt):]
        NFFT = 5000 
        Pxx, F = psd(sig - np.mean(sig), NFFT, Fs=1./dt,
                noverlap=int(NFFT*0.75), window=window_hanning)
        relP[nidx] = relativePower(Pxx, F, Frange=gRange)
        maxF[nidx] = maxPowerFrequency(Pxx, F, Frange=gRange)
        powerSpectra_P.append(Pxx)

    powerSpectra_F = F
    return relP, maxF, (np.array(powerSpectra_P), powerSpectra_F)



output_dir = 'output_local'
fileNamePrefix = ''
job_num = 0

in_fname = "{0}/{1}job{2:05}_output".format(output_dir, fileNamePrefix, job_num)
output_fname = in_fname
d = DataStorage.open(in_fname + ".h5", 'r')
options = d['options']
ei_net  = d['ei_net']

# Other
figSize = (12,8)
rcParams['font.size'] = 16
x_lim = [options['time']-2e3, options['time']]



###############################################################################


spikeMon_e = d['spikeMon_e']
spMon_ev_e   = spikeMon_e['status']['events']
senders_e    = spMon_ev_e['senders'] - spikeMon_e['gid_start']
spikeTimes_e = spMon_ev_e['times']

spikeMon_i = d['spikeMon_i']
spMon_ev_i = spikeMon_i['status']['events']
senders_i = spMon_ev_i['senders'] - spikeMon_i['gid_start']
spikeTimes_i = spMon_ev_i['times']

if ("pc_spikeMon" in d):
    PC = True
    pc_spikeMon = d['pc_spikeMon']
    pc_spMon_ev = pc_spikeMon['status']['events']
    senders_pc = pc_spMon_ev['senders'] - pc_spikeMon['gid_start']
    spikeTimes_pc = pc_spMon_ev['times']
else:
    PC = False

stateMonF_e = d['stateMonF_e']


events_e = []
events_i = []
events_e.append(d['stateMon_e'][0]['events'])
events_e.append(d['stateMon_e'][1]['events'])

events_i.append(d['stateMon_i'][0]['events'])
events_i.append(d['stateMon_i'][1]['events'])



F_tstart = 0.0
F_tend = options['time']
F_dt = 20.0
F_winLen = 250.0

Fe, Fe_t = slidingFiringRateTuple((senders_e, spikeTimes_e), ei_net['net_Ne'],
        F_tstart, F_tend, F_dt, F_winLen)

Fi, Fi_t = slidingFiringRateTuple((senders_i, spikeTimes_i), ei_net['net_Ni'],
        F_tstart, F_tend, F_dt, F_winLen)

Fpc, Fpc_t = slidingFiringRateTuple((senders_pc, spikeTimes_pc),
        ei_net['N_pc_created'], F_tstart, F_tend, F_dt, F_winLen)


bumpT = options['time'] - 2*F_winLen
bumpI = bumpT / F_dt
bump_e = np.reshape(Fe[:, bumpI], (ei_net['Ne_y'], ei_net['Ne_x']))
bump_i = np.reshape(Fi[:, bumpI], (ei_net['Ni_y'], ei_net['Ni_x']))

dim_e = Position2D()
dim_e.x = ei_net['Ne_x']
dim_e.y = ei_net['Ne_y']
G_est = fitGaussianBumpTT(bump_e, dim_e)
print G_est



# Flattened firing rate of E/I cells
figure(figsize=(12, 10))
subplot(3, 1, 1)
T, N_id = np.meshgrid(Fe_t, np.arange(ei_net['net_Ne']))
pcolormesh(T, N_id,  Fe)
#xlabel("Time (s)")
ylabel("Neuron #")
axis('tight')
colorbar()
title('Firing rate of E cells')
subplot(3, 1, 2)
T, N_id = np.meshgrid(Fi_t, np.arange(ei_net['net_Ni']))
pcolormesh(T, N_id,  Fi)
#xlabel("Time (s)")
ylabel("Neuron #")
axis('tight')
colorbar()
title('Firing rate of I cells')
# Flattened firing rate of Place cells
subplot(3, 1, 3)
T, N_id = np.meshgrid(Fpc_t, np.arange(ei_net['N_pc_created']))
xlim([0, 1000])
pcolormesh(T, N_id, Fpc)
xlabel("Time (s)")
ylabel("Neuron #")
axis('tight')
colorbar()
title('Firing rate of Place cells')
savefig(output_fname + '_FR_flat.png')


# External currents
figure()
ax = subplot(211)
plot(events_e[1]['times'], events_e[1]['I_stim'])
plot(events_e[0]['times'], events_e[0]['I_stim'])
ylabel('E cell $I_{stim}$')
axis('tight')
subplot(212)
plot(events_i[1]['times'], events_i[1]['I_stim'])
plot(events_i[0]['times'], events_i[0]['I_stim'])
ylabel('I cell $I_{stim}$')
xlabel('Time (ms)')


# E/I Vm
figure()
ax = subplot(211)
hold('on')
plot(events_e[1]['times'], events_e[1]['V_m'])
plot(events_e[0]['times'], events_e[0]['V_m'])
legend(['middle', 'edge'])
ylabel('E cell $V_m$')
subplot(212, sharex=ax)
plot(events_i[1]['times'], events_i[1]['V_m'])
plot(events_i[0]['times'], events_i[0]['V_m'])
legend(['middle', 'edge'])
ylabel('I cell $V_m$')
xlabel('Time (ms)')
xlim(x_lim)
savefig(output_fname + '_Vm.png')


# E/I I_syn
figure()
ax = subplot(211)
plot(events_e[1]['times'], events_e[1]['I_clamp_GABA_A'])
plot(events_e[0]['times'], events_e[0]['I_clamp_GABA_A'])
legend(['middle', 'edge'])
ylabel('E synaptic current (pA)')
subplot(212, sharex=ax)
plot(events_i[1]['times'], events_i[1]['I_clamp_AMPA'] + events_i[1]['I_clamp_NMDA'])
plot(events_i[0]['times'], events_i[0]['I_clamp_AMPA'] + events_i[0]['I_clamp_NMDA'])
legend(['middle', 'edge'])
xlabel('Time (s)')
ylabel('I synaptic current (pA)')
xlim(x_lim)
savefig(output_fname + '_Isyn.png')


# Firing rate of E cells on the twisted torus
figure()
pcolormesh(bump_e)
xlabel('E neuron no.')
ylabel('E neuron no.')
colorbar()
axis('equal')
title('Firing rates (torus) of E cells')
savefig(output_fname + '_firing_snapshot_e.png')


# Firing rate of I cells on the twisted torus
figure()
pcolormesh(bump_i)
xlabel('I neuron no.')
ylabel('I neuron no.')
colorbar()
axis('equal')
title('Firing rates (torus) of I cells')
savefig(output_fname + '_firing_snapshot_i.png')


### Firing rate of place cells on the twisted torus
#if (len(ei_net['PC']) != 0):
#    figure()
#    pcolormesh(np.reshape(Fpc[:, 0], (np.sqrt(ei_net['N_pc_created']),
#        np.sqrt(ei_net['N_pc_created']))))
#    xlabel('PC neuron no.')
#    ylabel('PC neuron no.')
#    colorbar()
#    axis('equal')
#    title('PC firing rates')


# Print a plot of bump position
(pos, bumpPos_times) = torusPopulationVector(
        (senders_e, spikeTimes_e), [ei_net['Ne_x'], ei_net['Ne_y']],
        tstart = options['theta_start_t'],
        tend   = options['time'],
        dt     = F_dt,
        winLen = F_winLen)
figure(figsize=figSize)
plot(bumpPos_times, pos)
xlabel('Time (s)')
ylabel('Bump position (neurons)')
legend(['X', 'Y'])
ylim([-ei_net['Ne_x']/2 -5, ei_net['Ne_y']/2 + 5])
savefig(output_fname + '_bump_position.pdf')


# Relative gamma power of inhibition in a sample of E neurons
Ftstart = 3 # sec
figure(figsize=(12, 6))
subplot(1, 2, 1)
gRange = (options['gammaRangeLow'], options['gammaRangeHigh'])
relP_e, maxF_e, spectra_e = relativeGammaPower(stateMonF_e, gRange,
        'I_clamp_GABA_A', tstart=Ftstart)
hist(relP_e)
xlabel('Relative power')
ylabel('Count')

subplot(1, 2, 2)
hist(maxF_e)
xlabel('Max. F (Hz)')
ylabel('Count')
savefig(output_fname + '_rel_gamma_P.pdf')


figure()
NP = 5
Fmax = 150
spectra_e_P = spectra_e[0]
spectra_e_F = spectra_e[1]
Frange = spectra_e_F <= Fmax
spectra_e_F = np.array([spectra_e_F]*spectra_e_P.shape[0])
plot(spectra_e_F[0:NP, Frange].T, spectra_e_P[0:NP, Frange].T)
xlabel('Frequency (Hz)')
ylabel('Power ($pA^2$/Hz)')
title('Power spectra of $I_{syn}$ of ' + str(NP) + ' selected E neurons')
savefig(output_fname + '_P_spectra.pdf')


#show()
