function [ipsc_id ipsc_sizes ipsc_start_id] = findIPSCs(sig, noise_std)
    % FINDEPSCS detect times and amplitudes of EPSCs from a voltage clamp
    % recording
    % Only EPSCs which have amplitude larger than 3*noise_std are
    % considered
   
    threshold = 4*noise_std;
    
    [xmax imax xmin imin] = extrema(sig);
    
    
    max_flag =  1 + zeros(1, numel(xmax));
    min_flag = -1 + zeros(1, numel(xmin));
    
    ex_ids = [imax imin];
    flags = [max_flag min_flag];
    
    [ex_ids flag_ids] = sort(ex_ids);
    flags = flags(flag_ids);
    
    ipsc_ex_ids = find(flags == 1);
    
    ipsc_id = ex_ids(ipsc_ex_ids);
    ipsc_start_id = ex_ids(find(flags == -1));
    ipsc_sizes = [0 (sig(ipsc_id(2:end)) - sig(ex_ids(ipsc_ex_ids(2:end) - 1)))];

    ipsc_thr = find(ipsc_sizes < threshold);
    ipsc_id(ipsc_thr) = [];
    ipsc_sizes(ipsc_thr) = [];
    
end