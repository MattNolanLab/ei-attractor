function r = lognrnd_local(mu,sigma, N1, N2)

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

r = exp(randn(N1, N2) .* sigma + mu);
