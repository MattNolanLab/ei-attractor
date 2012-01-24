function [vals] = sampleExponential(mean, N)
    % SAMPLEEXPONENTIAL generate N pseudorandom numbers drawn from an
    % exponential distribution with parameter 'mean'

    vals = -mean*log(1 - rand(1, N));
end