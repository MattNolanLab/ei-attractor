function [F, coherence] = getSigFreqCoherence(signal, dt, nPeriods)

    ddC_eps = 0;
    min_peak_ratio = 0;

    s_len = numel(signal);
    C = xcorr(signal - mean(signal), 'coeff');
    C = C(s_len:end);
    
    dC = diff(C);
    ddC = diff(dC);
    
    dC_zero = sign(dC(2:end) .* dC(1:end-1)) == -1;
    peak_i = find(dC_zero & ddC < ddC_eps) + 1;
    %peak_i(1) = [];
    C_peak = C(peak_i);
    
    % remove peaks which are insignificant compared to maximal peak
    [max_peak max_peak_i] = max(C_peak);
    
    if (numel(max_peak) == 0)
        F = 0
        coherence = 0;
    else
        t_max = (peak_i(max_peak_i)-1)*dt;
        F = 1/t_max;
        coherence = C((peak_i(max_peak_i)) * nPeriods);
    end
end