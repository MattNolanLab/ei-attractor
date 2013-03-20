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

from analysis.spikes    import slidingFiringRateTuple
from analysis.image     import Position2D, fitGaussianBumpTT
from analysis.signal    import relativePower, maxPowerFrequency
from data_storage       import DataStorage
from plotting.bumps     import bumpPosition, torusFiringRate, flatFiringRate, \
                               EIPlot



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
flatFiringRate(Fe, Fe_t, labelx="", titleStr = 'E cells')
subplot(3, 1, 2)
flatFiringRate(Fi, Fi_t, labelx="", titleStr = 'I cells')
# Flattened firing rate of Place cells
subplot(3, 1, 3)
flatFiringRate(Fpc, Fpc_t, labelx=None, titleStr = 'Place cells')
tight_layout()
savefig(output_fname + '_FR_flat.png')


# External currents
figure()
EIPlot(E = (events_e[1]['I_stim'], events_e[1]['times']),
       I = (events_i[1]['I_stim'], events_i[1]['times']),
       labely = "$V_m$")
EIPlot(E = (events_e[0]['I_stim'], events_e[0]['times']),
       I = (events_i[0]['I_stim'], events_i[0]['times']),
       labely = "$I_{stim}$")
xlim(x_lim)
savefig(output_fname + '_I_stim.png')


# E/I Vm
figure()
EIPlot(E = (events_e[1]['V_m'], events_e[1]['times']),
       I = (events_i[1]['V_m'], events_i[1]['times']),
       labely = "$V_m$")
EIPlot(E = (events_e[0]['V_m'], events_e[0]['times']),
       I = (events_i[0]['V_m'], events_i[0]['times']),
       labely = "$V_m$")
xlim(x_lim)
savefig(output_fname + '_Vm.png')


# E/I I_syn
figure()
EIPlot(E = (events_e[1]['I_clamp_GABA_A'], events_e[1]['times']),
       I = (events_i[1]['I_clamp_AMPA'] + events_i[1]['I_clamp_NMDA'], events_i[1]['times']),
       labely = "$I_{syn}$")
EIPlot(E = (events_e[0]['I_clamp_GABA_A'], events_e[0]['times']),
       I = (events_i[0]['I_clamp_AMPA'] + events_i[1]['I_clamp_NMDA'], events_i[0]['times']),
       labely = "$I_{syn}$")
xlim(x_lim)
savefig(output_fname + '_Isyn.png')


# Firing rate of E cells on the twisted torus
figure()
torusFiringRate(
        rateMap  = bump_e,
        labelx   = 'E neuron #',
        titleStr = 'Firing rates (torus) of E cells')
savefig(output_fname + '_firing_snapshot_e.png')


# Firing rate of I cells on the twisted torus
figure()
torusFiringRate(
        rateMap  = bump_i,
        labelx   = 'I neuron #',
        titleStr = 'Firing rates (torus) of I cells')
savefig(output_fname + '_firing_snapshot_i.png')


# Print a plot of bump position
figure()
bumpPosition(
    spikes    = (senders_e, spikeTimes_e),
    sheetSize = [ei_net['Ne_x'], ei_net['Ne_y']],
    tstart    = options['theta_start_t'],
    tend      = options['time'],
    dt        = F_dt,
    winLen    = F_winLen,
    units     = "ms")
savefig(output_fname + '_bump_position.pdf')



#show()
