% Script to transform spikeMonitor_times_n* to cell array

sheet_size = double(sheet_size);
clear spikeCell;

for nid = 0:sheet_size^2-1
        spikeCell{nid+1} = eval(['spikeMonitor_times_n' int2str(nid)]);
end