function [epsc_id epsc_sizes epsc_start_id] = findEPSCs(sig, noise_std)
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
    
    epsc_ex_ids = find(flags == -1);
    
    epsc_id = ex_ids(epsc_ex_ids);
    epsc_start_id = ex_ids(find(flags == 1));
    epsc_sizes = -(sig(epsc_id) - sig(ex_ids(epsc_ex_ids - 1)));

    epsc_thr = find(epsc_sizes < threshold);
    epsc_id(epsc_thr) = [];
    epsc_sizes(epsc_thr) = [];
    
end