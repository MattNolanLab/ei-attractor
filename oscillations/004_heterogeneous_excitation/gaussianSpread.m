function [Y] = gaussianSpread(X, sigma, maxVal)
    Y = maxVal*exp(-X.^2/sigma^2);
end