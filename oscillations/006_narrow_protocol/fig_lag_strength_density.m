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

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);


mean_first_e = zeros(1, nParam);
mean_first_i = zeros(1, nParam);

for par_it = 1:nParam
    res = results(par_it, trial_it);
    D = res.opt.D;
    dt = res.opt.dt;

    first_spike_e = zeros(1, Ne);
    first_spike_i = zeros(1, Ni);

    for n_it = 1:Ne
        if numel(res.spikeCell_e{n_it}) ~= 0
            first_spike_e(n_it) = res.spikeCell_e{n_it}(1);
        end
    end
    first_spike_e(find(first_spike_e == 0)) = [];
    mean_first_e(par_it) = mean(first_spike_e);


    for n_it = 1:Ni
        if (numel(res.spikeCell_i{n_it}) ~= 0)
            first_spike_i(n_it) = res.spikeCell_i{n_it}(1);
        end
    end        
    first_spike_i(find(first_spike_i == 0)) = [];
    mean_first_i(par_it) = mean(first_spike_i);
end

we_vec = res.opt.we_vec;
sp_vec = res.opt.sparseness_vec;

mean_first_e = reshape(mean_first_e, Nwe, Nsp);
%std_first_e =  reshape(std(first_spike_e'), Nwe, Nsp);

mean_first_i = reshape(mean_first_i, Nwe, Nsp);
%std_first_i = reshape(std(first_spike_i'), Nwe, Nsp);


figure;
[grid_sp, grid_we] = meshgrid(sp_vec, we_vec);
surf(grid_sp, grid_we, mean_first_i - mean_first_e);
xlabel('Sparseness');
ylabel('Synaptic strength');
