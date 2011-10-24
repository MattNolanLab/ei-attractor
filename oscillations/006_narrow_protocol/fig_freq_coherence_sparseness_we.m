% Plot frequency of oscillation and coherence as a function of sparseness
% vs. syn. strength

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '010';

nParam  = size(results, 1);
nTrials = size(results, 2);

mA = 1000;
fontSize = 14;

spread_all = [0.5 2 5 10];

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

sp_vec = results(1,1).opt.sparseness_vec;
we_vec = results(1,1).opt.we_vec;

par_it = 20;
trial_it = 1;

Ni_it = 1;
Ne_it = 1;


F = 10:1:100;
sampling_rate = 1e4;
win_len = 2000;
noverlap = win_len /2;
plot_dB = false;

dt = results(1,1).opt.dt;
nPeriods = 1;
t_start = 1;
t_start_i = fix(t_start/dt) + 1;


freq = zeros(Nsp, Nwe);
coherence  = zeros(Nsp, Nwe);
%freq_std  = zeros(Nsp, Nwe);
%pow_mean = zeros(Nsp, Nwe);

par_it = 1;
for sp_it = 1:Nsp
    for we_it = 1:Nwe

        res = results(par_it, trial_it);
        

        ge_sig = res.Vmon.ge(Ni_it, t_start_i:end)*mA;
        gi_sig = -res.Vmon.gi(Ne_it, t_start_i:end)*mA;

%         % Stellate cell
%         [Y, F, T, P] = spectrogram(gi_sig, win_len, noverlap, F, sampling_rate);
%         % The following code produces the same result as calling
%         % spectrogram with no outputs:
% 
% 
%         if (plot_dB)
%             P_plot = 10*log10(abs(P));
%             power_label = 'Power (dB)';
%         else
%             P_plot = abs(P);
%             power_label = 'Power';
%         end
% 
%         [max_P max_P_i] = max(P_plot);
%         fmax = F(max_P_i);
% 
%         freq_mean(sp_it, we_it) = mean(fmax);
%         pow_mean(sp_it, we_it) = mean(max_P);
%         freq_std(sp_it, we_it) = std(fmax);

        [F coh] = getSigFreqCoherence(gi_sig, dt, nPeriods);
        
        freq(sp_it, we_it) = F;
        coherence(sp_it, we_it) = coh;
        
        par_it = par_it+1;
    end
end

figure('Position', [800 1050 900 1050]);
[WE SP] = meshgrid(we_vec*1000, sp_vec);
subplot(2,1,1);
surf(SP, WE, freq);
ylabel('Syn. strength (mA)');
xlabel('Sparseness');
view(-64, 46);

subplot(2,1,2);
surf(SP, WE, coherence);

ylabel('Syn. strength (mA)');
xlabel('Sparseness');
view(-64, 46);


set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print('-depsc2', sprintf('%s/%s_freq_N%.3d.eps', ...
        outputDir, outputNum, Ne_it));    



