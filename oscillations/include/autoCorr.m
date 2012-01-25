function C = autoCorr(signal)
    s_len = numel(signal);
    C = xcorr(signal - mean(signal), 'coeff');
    C = C(s_len:end);
end