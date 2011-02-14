% Print coloured 2d histogram of the population response
% function printPopulationResponse(sheet_size, startTime, endTime)

%close all;

    opt = parseOptions(options);
    sheet_size = opt.sheet_size;

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    delta_t = 1; % sec
    startTime = 17;
    endTime = startTime; % sec
    
    firingPop = zeros(sheet_size, sheet_size);
    
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;
            neuronSpikes = spikeCell{neuronID + 1}; %eval(['spikeMonitor_times_n' num2str(neuronID)]);
            firingRate = computeFiringRate(neuronSpikes, startTime, endTime, dt_rat, delta_t);
            
            firingPop(x_i+1, y_i+1) = firingRate(numel(firingRate));
        end
    end
    
    firingPop = firingPop';

    %histmat = hist2(id_x, id_y, xedges, yedges);
    figure(3);
    hold off;
    pcolor(0:sheet_size-1,0:sheet_size-1,firingPop);
    axis square tight;
    colorbar;
    shading flat;
    colormap(jet);
    
    %SNList_nID = 38;
    %drawPin(double(SNList(SNList_nID)), sheet_size, [1 1 0]);
    
%     fourierPop = rot90(fftshift(fft2(firingPop)));
%     absFourierPop = abs(fourierPop);
%     figure(2);
%     pcolor(0:sheet_size-1, 0:sheet_size-1, absFourierPop);
%     
%     colormap(gray);
%     colorbar;
%     axis square tight;
%     %shading interp;
%     
%     [maxF maxF_y] = max(absFourierPop);
%     [sortedF sorted_I] = sort(maxF, 'descend');
%     
%     maxF1 = sortedF(1);
%     maxF1_x = sorted_I(1);
%     maxF1_y = maxF_y(maxF1_x);
%     
%     maxF2 = sortedF(2);
%     maxF2_x = sorted_I(2);
%     maxF2_y = maxF_y(maxF2_x);
%    
%     max1Vec = complex(maxF1_x, maxF1_y);
%     max2Vec = complex(maxF2_x, maxF2_y);
%     an = angle(max1Vec-max2Vec)/2/pi * 360
%     
%     normF = absFourierPop/max(max(absFourierPop));
%     thrF = double(im2bw(normF, 1/2));
%     figure(3);
%     pcolor(0:sheet_size-1, 0:sheet_size-1, thrF);
%  
%     colormap(gray);
%     colorbar;
%     axis square tight;
%     %shading interp;
%     
%     
%     figure(4);
%     firingThr = 0.35;
%     thrFiringPop = double(im2bw(firingPop/max(max(firingPop)), firingThr))
%     pcolor(thrFiringPop);
%     axis square tight;
%     colorbar;
%end
