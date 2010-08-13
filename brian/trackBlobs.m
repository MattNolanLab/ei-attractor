function [segCenters_r segCenters_c, thrFiringPop, segFiringPop] = trackBlobs(firingPop)

    sheet_size = size(firingPop, 1);
    
    % Simply threshold the population response to segment the image
    % This should easily work, since the blobs are coherent
    firingThr = 0.35;
    blobNThreshold = 15;
    
    thrFiringPop = zeros(sheet_size);
    thr_i = find(firingPop/max(max(firingPop)) >= firingThr);
    %thrFiringPop = reshape(firintPop, sheet_size^2, 1);
    thrFiringPop(thr_i) = 1;
    thrFiringPop = reshape(thrFiringPop, sheet_size, sheet_size);

    % Segment no. for each pixel
    % -1: no segment
    % >0: segment number
    segFiringPop = -1 + zeros(sheet_size);
    
    % Find nonzero elements, i.e. blobs in the matrix and identify clusters
    [nz_r, nz_c] = find(thrFiringPop);
    
    segn = 0; % segment no. iterator    
    for it = 1:numel(nz_r)
        
        if (segFiringPop(nz_r(it), nz_c(it)) == -1)  % new blob
            visitPoints_r = [nz_r(it)];
            visitPoints_c = [nz_c(it)];
            
            while numel(visitPoints_r) ~= 0
                % find neighbours, 8-neighborhood
                r = visitPoints_r(1); visitPoints_r(1) = [];
                c = visitPoints_c(1); visitPoints_c(1) = [];
                segFiringPop(r,c) = segn;
                
                neighb_r = [r-1, r-1, r-1, r,   r,   r+1, r+1, r+1];
                neighb_c = [c-1, c,   c+1, c-1, c+1, c-1, c,   c+1];
                
                for n_i = 1:8
                    if (neighb_r(n_i) >= 1 && neighb_r(n_i) <= sheet_size && ...
                        neighb_c(n_i) >= 1 && neighb_c(n_i) <= sheet_size && ...
                        thrFiringPop(neighb_r(n_i), neighb_c(n_i)) == 1 && ...
                        segFiringPop(neighb_r(n_i), neighb_c(n_i)) == -1)
                        
                        visitPoints_r = [visitPoints_r neighb_r(n_i)];
                        visitPoints_c = [visitPoints_c neighb_c(n_i)];
                        segFiringPop(neighb_r(n_i), neighb_c(n_i)) = segn;
                    end
                end
            end
            
            segn = segn+1;
        end
    end

    % Search through the clusters and find the centers of blobs
    segCenters_r = [];
    segCenters_c = [];
    for it = 1:segn
        [r, c] = find(segFiringPop == it-1);
        
        % reject blobs which have less than the threshold number of points
        if (numel(r) >= blobNThreshold)
            segCenters_r = [segCenters_r mean(r)];
            segCenters_c = [segCenters_c mean(c)];
        end
    end
end
