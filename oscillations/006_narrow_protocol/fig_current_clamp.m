% Plot lag of onset of spiking between interneurons and stellate cells, as
% a function of excitatory synaptic strength x density (sparsity)

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '008';

nParam  = size(results, 1);
nTrials = size(results, 2);

mA = 1000;
fontSize = 14;

spread_all = [0.5 2 5 10];

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

par_it = 1;
trial_it = 1;

Ne_it = 10;
Ni_it = 1;

res = results(par_it, trial_it);

% % E and I synaptic kernels
% ker_dt = res.opt.dt;
% 
% I_ker_T = 10 * res.opt.taui;
% I_ker_t = 0:ker_dt:I_ker_T;
% I_ker = exp(-I_ker_t/res.opt.taui);
% 
% E_ker_T = 10 * res.opt.taue;
% E_ker_t = 0:ker_dt:E_ker_T;
% E_ker = exp(-E_ker_t/res.opt.taue);
% 
% % Plot synaptic currents onto stellate cell
% pre_I_cells = find(results(par_it, trial_it).Mi(Ne_it, :) == 1);
% if (numel(pre_I_cells) == 0)
%     error(sprintf('No I cells were connected to E cell no. %d', Ne_it));
% end


% pre_E_cells = find(res.Me(Ni_it, :) == 1);
% if (numel(pre_E_cells) == 0)
%     error(sprintf('No E cells were connected to I cell no. %d', Ni_it));
% end
% 
% figure;
% I_post_spikes = full(sum(res.spikeRecord_e(pre_E_cells, :))) * res.opt.we;
% I_post_current = conv(I_post_spikes, E_ker, 'same');
figure();
plot(res.times, res.Vmon.ge*1000);
xlabel('Time (s)');
ylabel('Syn. current (mA)');
title('Interneuron');



figure();

gi_sig = - res.Vmon.gi*mA;

% E_post_spikes = full(sum(res.spikeRecord_i(pre_I_cells, :))) * res.opt.wi;
% E_post_current = conv(E_post_spikes, I_ker, 'same');
subplot(2,1,1, 'FontSize', fontSize);
plot(res.times, gi_sig);
xlabel('Time');
ylabel('Syn. current (mA)');
title('Stellate cell');


set(gca, 'FontSize', fontSize);
F = 0:2:150;
[Y, F, T, P] = spectrogram(gi_sig, 1024, 940, F, 1/res.opt.dt);
% The following code produces the same result as calling
% spectrogram with no outputs:

%filt = find(P < 1e-1);
%P(filt) = nan;
plot_dB = false;

if (plot_dB)
    P_plot = 10*log10(abs(P));
    power_label = 'Power (dB)';
else
    P_plot = abs(P);
    power_label = 'Power';
end

min_P = min(min(P_plot));
max_P = max(max(P_plot));

subplot(2,1,2, 'FontSize', fontSize)
surf(T,F,P_plot,'EdgeColor','none');
axis xy; axis tight; %view(0,90);
%zlim([min_P max_P]);
xlabel('Time');
ylabel('Frequency (Hz)');
zlabel(power_label);
colormap jet;
view(-90, -90)
%colorbar;

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'wide_c003_Current_2_spectrogram_non_dB.eps');


