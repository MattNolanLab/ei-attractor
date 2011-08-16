function [cc, edges] = crossCorrelation(N1, N2, dt, range, T)
% CROSSCORRELATION Compute correlations statistics between two neurons

cc = ISIpairs(N1, N2);

edges = [fliplr((dt/2:dt:range)*-1) dt/2:dt:range];
cc = histc(cc, edges);
%cc = cc - sz*dt/T^2;
bar(gca, edges, cc, 'histc')
grid on;
axis tight;

xlabel('Time (s)');
ylabel('Interval count');

end