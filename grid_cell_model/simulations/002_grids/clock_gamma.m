%
%   clock_gamma.m
%
%   Gamma as a clock script.
%
%       Copyright (C) 2012  Hugh Pastoll
%       
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
%       
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
%       
%       You should have received a copy of the GNU General Public License
%       along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
close all

fontSize = 15


fs = 10000; % sampling frequency
dt = 1./fs;
bumpPos_dt = bumpPos_times(2) - bumpPos_times(1)

thetafreq = 8; % light stimulation frequency

N = size(theta_stateMon_Iclamp_e_values, 1)

trace = theta_stateMon_Iclamp_e_values(2, :); % the data


stim_dur = .5; % stimulus duration
stimstart = 0; % time of stimulus onset in recording

pa = 1e12; % scaling factor to convert amps to picoamps

% Filtering

hd = bandpassgamma(fs); % needs a helper function
sfilt = filtfilthd(hd, trace*pa); % in pA

sfilt = sfilt - mean(sfilt); % zero mean

% extract the stimulation
trace = sfilt(stimstart*fs+1:(stimstart+stim_dur)*fs)';
overlay = reshape(trace, fs/thetafreq, []);

% remove first theta cycle to eliminate stimulus onset artifacts
%overlay(:, 1) = [];
overlay_mean = mean(overlay, 2); % average episode

% use variance of first and last quarter of each episode as baseline
pts = round(fs/thetafreq/8); % number of point after beginning and before end

% find the initial maxima of the trace that exceeds 3 std devs of the first
% and last portions
[pks, locs] = findpeaks(overlay_mean, 'minpeakheight', 3*std(overlay_mean([1:pts end-pts:end])));
t = linspace(0, 1/thetafreq, length(overlay_mean));
trigger = t(locs(1)); % phase drift is relative to this time point


% use entire theta cycle (as opposed to some portion)
active_gamma = overlay;

% calculate phase
ht = hilbert(active_gamma); % hilbert transform 
phase = angle(ht); % extract phase information from the hilbert transform
drift = unwrap(phase); % cumulative phase

% matrix for storing the sorted cumulative phase (phase drift)
drift_sorted = zeros(size(drift)); 

% point at which to collect phases into the same phase bracket
ref_point = locs(1); % first gamma peak > 3 std devs above the mean
ref = 0;

% Line up unwrapped phases at a point by shifting up or down by a
% factor of 2 pi
const = 2*pi;

% interate through all traces
for i = 1:size(overlay, 2)
    
    fac = 0;
    val = drift(ref_point, i);
    % while there is still a bigger than pi difference between the phase mean and
    % the phase value adjust the phase factor in the direction of the difference
    while abs(ref-(val+fac*const))>const/2
        
        fac = fac + sign(ref-val);
    end
    
    % add the conversion factor to the original unwrapped phase matrix
    drift_sorted(:, i) = drift(:, i) + fac*const;
    
end

spread = std(drift_sorted, [], 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bump position versus gamma phase at reference point
Ntheta = size(overlay, 1);
Nepochs = size(overlay, 2);
gamma_ref_t = stimstart + (ref_point + Ntheta * [0:Nepochs-1])*dt;
bumpPos_at_ref = fix(gamma_ref_t/bumpPos_dt)

% actual bump position
d_bump = sqrt(bumpPos(bumpPos_at_ref, 1).^2 + bumpPos(bumpPos_at_ref, 2).^2);
%d_bump = bumpPos(bumpPos_at_ref, 1)
figure
subplot(1,1,1, 'FontSize', fontSize)
plot(d_bump, drift_sorted(ref_point, :), 'o')
xlabel('Bump distance from center (neurons)')
ylabel('\gamma-phase at ref. point')

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'output_local/clock_gamma_phase_bump_pos.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


example_it = round(size(overlay, 2)/2)

% %% plotting
figure('Position', [0, 0, 600 1050])
m = 5; n = 1;
subplot(m, n, 1, 'FontSize', fontSize)
p = pcolor(overlay');
set(p,'EdgeColor','none')
set(gca, 'XTick', [])
ylabel('\theta epoch')

subplot(m, n, 2, 'FontSize', fontSize)
plot(overlay(:, example_it))
ylabel('Current (pA)')
xlim([0 size(overlay, 1)])
set(gca, 'XTick', [])

subplot(m, n, 3, 'FontSize', fontSize)
plot(phase(:, example_it))
ylabel('\gamma-phase')
xlim([0 size(overlay, 1)])
set(gca, 'XTick', [])

subplot(m, n, 4, 'FontSize', fontSize)
plot(drift_sorted, 'Color', [0.5 0.5 0.5])
hold on
plot(drift_sorted(:, example_it), 'k', 'Linewidth', 2)
hold off
ylabel('cum. \gamma-phase')
set(gca, 'XTick', [])
xlim([0 size(overlay, 1)])

subplot(m, n, 5, 'FontSize', fontSize)
plot(spread, 'k--', 'LineWidth', 1.5)
ylabel('\sigma  cum. \gamma-phase')
set(gca, 'XTick', [0 size(overlay, 1)])
set(gca, 'XTickLabel', {'-3.14', '3.14' })
xlim([0 size(overlay, 1)])

ph_noise = std(2*pi*rand(size(overlay, 2), 1))
hold on
plot(0:size(overlay, 1)-1, ph_noise+zeros(size(overlay, 1), 1), '--r')
ylim([0 4])

set(gcf,'PaperPositionMode','auto');
print('-djpeg90', '-r300', 'output_local/clock_gamma.jpg')


figure('Position', [0, 0, 800 400])
m = 1; n = 2;
subplot(m, n, 1, 'FontSize', fontSize)
p = pcolor(overlay');
set(p,'EdgeColor','none')
set(gca, 'XTick', [])
ylabel('\theta epoch')
xlabel('\theta phase')
set(gca, 'XTick', [1 size(overlay, 1)])
set(gca, 'XTickLabel', {'-3.14', '3.14' })

subplot(m, n, 2, 'FontSize', fontSize)
d_bump_all = sqrt(diff(bumpPos(:, 1))/bumpPos_dt.^2 + diff(bumpPos(:, 2))/bumpPos_dt.^2);
%d_bump_all = bumpPos(:, 1)

plot(d_bump_all, [0:numel(bumpPos_times)-2]*bumpPos_dt)
ylim([0 stim_dur])
%xlabel('Bump distance from center of the sheet (neurons)')
xlabel('Bump velocity (neur./s)')
ylabel('Time (s)')

set(gcf,'PaperPositionMode','auto');
print('-djpeg80', '-r300', 'output_local/clock_theta_bump_pos.jpg')

