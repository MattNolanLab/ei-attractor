function [F] = getPopOscillationFreq(firingRate, dc_ratio, dt)
    % Find peak of oscillation in the firing rate signal. If ratio
    %  f_peak/f_0 < dc_ratio, the resulting F will be 0 (asynchronous mode)
    %  f_peak is computed from signal with DC component subtracted
    %  f_0 is from signal with DC component
    
    [Y f NFFT] = fourierTrans(firingRate, dt);
    Y_abs = 2*abs(Y(1:NFFT/2+1));

    sig_0 = firingRate - mean(firingRate);
    [Y_0 f_0 NFFT_0] = fourierTrans(sig_0, dt);
    Y_0_abs = 2*abs(Y_0(1:NFFT_0/2+1));
    

    [maxF maxFI] = max(Y_abs);
    [maxF_0 maxFI_0] = max(Y_0_abs);
    
    if (maxF_0/maxF < dc_ratio)
        F = nan;
    else
        F = f_0(maxFI_0);
    end

end