function [isi] = ISIpairs(N1, N2)
% ISIPAIRS create a vector os ISI pairs of neurons N1 and N2, containing
% times of spikes
%
% Absolute value not provided

[mesh1, mesh2] = meshgrid(N1, N2);
isi = mesh1 - mesh2;
isi = reshape(isi, 1, numel(N1)*numel(N2));