function [F, coherence] = getPopOscFreqAutoCorr(signal, dt, nPeaks)

    ddC_eps = 0;
    min_peak_ratio = 1/10;

    s_len = numel(signal);
    C = xcorr(signal - mean(signal), 'coeff');
    C = C(s_len:end);
    
    dC = diff(C);
    ddC = diff(dC);
    
    dC_zero = sign(dC(2:end) .* dC(1:end-1)) == -1;
    peak_i = find(dC_zero & ddC < ddC_eps);
    C_peak = C(peak_i);
    
    % remove peaks which are insignificant compared to maximal peak
    max_peak = max(C_peak);
    peak_i(find(C_peak < min_peak_ratio)) = [];
    
    if numel(peak_i) >= nPeaks
        % Add peak at zero to filtered peaks
        peak_i = [1 peak_i];
        C_peak = C(peak_i);

        T = mean(peak_i(2:nPeaks) - peak_i(1:nPeaks-1));
        F = 1 / T / dt;
        coherence = mean(C_peak(2:nPeaks));
    else
        % no peaks detected, assume asynchronous signal
        F = 0;
        coherence = 0;
    end
end