function [G crossCorr angles] =  cellGridness(firingHist, edges)

cutR = 90; %cm
ac_cutRmin = 20;
autoCorr_edges = linspace(edges(1)*2, edges(end)*2, numel(edges)*2 - 1);

firingHist(isnan(firingHist)) = 0;
autoCorr = xcorr2(firingHist - mean(reshape(firingHist, 1, numel(firingHist))));

[X Y] = meshgrid(autoCorr_edges);
autoCorr(sqrt(X.^2 + Y.^2) < ac_cutRmin) = 0;

%autoCorr(sqrt(X.^2 + Y.^2) > cutR) = 0;
% autoCorr(find(abs(autoCorr_edges) > cutR), :) = [];
% autoCorr(:, find(abs(autoCorr_edges) > cutR)) = [];
% autoCorr = autoCorr - mean(reshape(autoCorr, 1, numel(autoCorr)));

it=1;
angles = 0:3:180;
crossCorr = zeros(1, numel(angles));
for angle = angles
    autoCorrRot = imrotate(autoCorr, -angle, 'bicubic', 'crop');
    C = corrcoef(reshape(autoCorr, 1, numel(autoCorr)), ...
        reshape(autoCorrRot, 1, numel(autoCorr)));
    crossCorr(it) = C(1, 2);
    
    it = it+1;
end

maxima = max(crossCorr([11 31 51])); %30, 90, and 150 deg
minima = min(crossCorr([21 41])); % 60 and 120 deg

G = minima - maxima;

end