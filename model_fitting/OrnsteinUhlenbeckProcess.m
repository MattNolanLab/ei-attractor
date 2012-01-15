function [X] = OrnsteinUhlenbeckProcess(tau, mju, sigma, dt, T)
% Ornstein-Uhlenbeck process

t = 0:dt:T;
X = zeros(1,numel(t));

X(1) = mju;
for it = 1:numel(t)-1
    X(it+1) = X(it) + (mju - X(it))/tau*dt + sqrt(2*sigma^2*dt/tau)*randn();
end


end
