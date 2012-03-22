function hd = bandpassgamma(fs)
% Bandpassintra returns a discrete-time filter object.


% parameters (cf. passband cutoffs from Chrobak & Buzsaki 1998 50-150, with a very
% approx 120dB/octave roll off)

% All frequency values are in Hz.
Fs = fs;  % Sampling Frequency

% NB Do not make parameters more extreme - they will make the filter
% unstable

% was 10, 15 fstop1 and fpass1

Fstop1 = 40;     % First Stopband Frequency
Fpass1 = 50;     % First Passband Frequency
Fpass2 = 120;  % Second Passband Frequency
Fstop2 = 150;  % Second Stopband Frequency
Astop1 = 10;     % First Stopband Attenuation (dB)
Apass  = 1;     % Passband Ripple (dB)
Astop2 = 10;     % Second Stopband Attenuation (dB)

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord([Fpass1 Fpass2]/(Fs/2), [Fstop1 Fstop2]/(Fs/2), Apass, ...
                 max(Astop1, Astop2));

             


             
% Calculate the zpk values using the BUTTER function.
[z,p,k] = butter(N, Fc);

% To avoid round-off errors, do not use the transfer function.  Instead
% get the zpk representation and convert it to second-order sections.
[sos_var,g] = zp2sos(z, p, k);
hd          = dfilt.df2sos(sos_var, g);

% [EOF]
