function [e_mfr_all, i_mfr_all] = meanFiringRateAll(res, t_start, t_end)

        t_start_i = t_start/res.opt.dt + 1;
        t_end_i   = t_end/res.opt.dt + 1;
    
        times = res.times(t_start_i:t_end_i);

        % mean firing rates of neurons in this trial
        mfr_T = t_end - t_start;
        
        e_mfr_all = full(sum(res.spikeRecord_e')/mfr_T);
        i_mfr_all = full(sum(res.spikeRecord_i')/mfr_T);
end