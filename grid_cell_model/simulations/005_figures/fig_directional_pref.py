#
#   fig_directional_pref.py
#
#   Measure of directional preference of model grid cells (with/without place
#   cell input)
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

from scipy.io           import loadmat
from matplotlib.pyplot  import *
from tables             import *

from analysis.grid_cells import SNFiringRate, motionDirection


jobNum = 3100
trialNum = 0
dumpNum = 7


rcParams['font.size'] = 14

arenaDiam = 180.0     # cm
h = 3.0
Ne = 34

# Neuron to extract spikes from
neuronNum = 1
spikeType = 'excitatory'
winLen = 0.5 # s
tend = 1200 #s
binWidth = 6.0/360.0 * 2*np.pi

directionalN_start = 68   # Total number of resulting directional vectors plotted
directionalN_end   = 87
directionalRange = range(0, 20) + range(68, 88)


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"


def binRates(rates, angles, binWidth, dt):
    '''
    Bin firing rates into bins of angles of width binWidth (radians).
    dt is the time sampling - used for counting relative time spent in each bin
    '''
    angles_shifted = angles + np.pi
    bin_angles  = np.arange(0, 2*np.pi, binWidth)
    rates_mean  = np.ndarray(len(bin_angles))
    rates_count = np.ndarray(len(bin_angles))
    bin_its     = np.floor(angles_shifted/binWidth)
    for it in xrange(len(bin_angles)):
        rates_mean[it] = np.mean(rates[bin_its == it])
        rates_count[it] = np.sum(bin_its == it) * dt
    return rates_mean, bin_angles, rates_count


def meanDirectionVector(angles, rates):
    '''
    Simply calculate a mean vector, taken from polar coordinates
    Return in polar coordinates.
    '''
    x = np.mean(np.cos(angles) * rates)
    y = np.mean(np.sin(angles) * rates)
    return np.arctan2(y, x), np.sqrt(x**2 + y**2)


def joinPolarPlotData(data):
    return np.append(data, data[0])



fileName = fileNameTemp
fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                            trialNum, dumpNum)
try:
    data = loadmat(fileName +  '_output.mat')
except:
    print "warning: could not open: " + fileName
    exit(1)

pos_x           = data['pos_x'].ravel()
pos_y           = data['pos_y'].ravel()
rat_dt          = data['dt'][0][0]
velocityStart   = data['velocityStart'][0][0]
if spikeType == 'excitatory':
    spikeTimes  = data['spikeCell_e'].ravel()
if spikeType == 'inhibitory':
    spikeTimes  = data['spikeCell_i'].ravel()

gridSep         = data['options']['gridSep'][0][0][0][0]


spikes = spikeTimes[neuronNum] - velocityStart*1e-3
spikes = np.delete(spikes, np.nonzero(spikes < 0)[0])

rate, rate_t = SNFiringRate(spikes, tend, rat_dt, winLen)
angles, angle_t, avg_spd = motionDirection(pos_x, pos_y, rat_dt, tend, winLen)

figure()
plot(rate_t, rate)
figure()
plot(angle_t, angles)

figure()
polar(angles, rate, 'o')

rates_mean, bin_angles, rates_count = binRates(rate, angles, binWidth, rat_dt)
bin_angles_join = joinPolarPlotData(bin_angles)
rates_mean_join = joinPolarPlotData(rates_mean)
rates_count_join = joinPolarPlotData(rates_count)
md_angle, md_r = meanDirectionVector(bin_angles, rates_mean)


figure(figsize=(6, 10))
subplot(2, 1, 1, projection="polar")
polar(bin_angles_join, rates_mean_join)
annotate("", xytext=(0.0,0.0), xy=(md_angle, md_r), arrowprops=dict(facecolor='black',
    width=1, headwidth=6, frac=0.2))
title('Average firing rate per bin (Hz)')
tight_layout()
subplot(2, 1, 2, projection="polar")
polar(bin_angles_join, rates_count_join)
title('Time spent in each bin (s)')
tight_layout()

savefig(fileName + '_dir_tuning.pdf')


figure()
subplot(1, 1, 1, polar=True)
# Now do all the things for some of the neurons in the population (e.g. first 20)
# and only plot the resultant average directional vector
for n_it in directionalRange:
    print("Vector no. " + str(n_it))
    spikes = spikeTimes[n_it] - velocityStart*1e-3
    spikes = np.delete(spikes, np.nonzero(spikes < 0)[0])
    
    rate, rate_t = SNFiringRate(spikes, tend, rat_dt, winLen)
    angles, angle_t, avg_spd = motionDirection(pos_x, pos_y, rat_dt, tend, winLen)

    rates_mean, bin_angles, rates_count = binRates(rate, angles, binWidth, rat_dt)
    md_angle, md_r = meanDirectionVector(bin_angles, rates_mean)

    annotate("", xytext=(0.0,0.0), xy=(md_angle, md_r), arrowprops=dict(facecolor='black',
        width=1, headwidth=6, frac=0.2))

    ylim([0, 0.5])

savefig(fileName + '_dir_tuning_cells.pdf')

show()
