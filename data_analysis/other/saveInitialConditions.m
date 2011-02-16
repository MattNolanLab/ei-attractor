% Extract membane potential of all neurons from SNMonitor_values (at the
% end of the simulation) and save them with other definining parameters
% (e.g. options)

opt = parseOptions(options);

fileName = ['initial_conditions_' int2str(opt.sheet_size) '.mat'];

time_id = numel(SNMonitor_times);

snapshot_time = SNMonitor_times(time_id);
membrane_potentials = SNMonitor_values(:, time_id)';
synaptic_cond = SNgMonitor_values(:, time_id);

display(['Exporting initial conditions to file...' fileName]);
save(fileName, 'options', 'opt', 'snapshot_time', 'membrane_potentials', ...
    'synaptic_cond');