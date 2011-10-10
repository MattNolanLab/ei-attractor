% Plot lag of onset of spiking between interneurons and stellate cells, as
% a function of excitatory synaptic strength x density (sparsity)

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '005';

nParam  = size(results, 1);
nTrials = size(results, 2);

spread_all = [0.5 2 5 10];

%par_it = 1;
trial_it = 1;

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);


first_spike_e = zeros(nParam, Ne);
first_spike_i = zeros(nParam, Ni);

for par_it = 1:nParam
    res = results(par_it, trial_it);
    D = res.opt.D;
    dt = res.opt.dt;

    for n_it = 1:Ne
        if numel(res.spikeCell_e{n_it}) ~= 0
            first_spike_e(par_it, n_it) = res.spikeCell_e{n_it}(1);
        end
    end

    for n_it = 1:Ni
        if (numel(res.spikeCell_i{n_it}) ~= 0)
            first_spike_i(par_it, n_it) = res.spikeCell_i{n_it}(1);
        end
    end

end

we_vec = res.opt.we_vec;

mean_first_e = mean(first_spike_e');
std_first_e = std(first_spike_e');

mean_first_i = mean(first_spike_i');
std_first_i = std(first_spike_i');

figure;
bar(we_vec, mean_first_i - mean_first_e);
