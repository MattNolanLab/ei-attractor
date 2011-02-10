function [Y, f, NFFT] = fourierTrans(signal, dt)
% FOURIERTRANS Compute fast Fourier transform of signal, with time step dt

Fs = 1/dt;
signalMean = mean(signal);
signal = signal - signalMean;
sL = numel(signal);

NFFT = 2^nextpow2(sL); % Next power of 2 from length of y
Y = fft(signal,NFFT)/sL;
f = Fs/2*linspace(0,1,NFFT/2+1);

end