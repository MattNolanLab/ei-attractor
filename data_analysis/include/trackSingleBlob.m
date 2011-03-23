function [segCenters_r segCenters_c, thrFiringPop, segFiringPop] = trackSingleBlob(firingPop)
% TRACKSINGLEBLOB tracks a single blob on square sheet. Boundaries are
% wrapped around (torus)
%
% firingPop is an array (1D) of firing rates of the population. This will
% be reshaped to a square
%
% Here we are assuming that the population contains only one bump on a 2d
% torus

    sheet_size = sqrt(numel(firingPop));
    firingPop = reshape(firingPop, sheet_size, sheet_size)'; % reshape takes all by columns
    
    % Simply threshold the population response to segment the image
    % This should easily work, since the blobs are coherent
    firingThr = 0.35;
    blobNThreshold = 15;
    
    thrFiringPop = zeros(sheet_size);
    thr_i = find(firingPop/max(max(firingPop)) >= firingThr);
    thrFiringPop(thr_i) = 1;
    %thrFiringPop = reshape(thrFiringPop, sheet_size, sheet_size);

    % Segment no. for each pixel
    % -1: no segment
    % >0: segment number
    segFiringPop = -1 + zeros(3*sheet_size);
    
    % Find nonzero elements, i.e. blobs in the matrix and identify clusters
    [nz_r, nz_c] = find(thrFiringPop);
    nz_r = nz_r + sheet_size;
    nz_c = nz_c + sheet_size;
    
    segn = 0; % segment no. iterator    
    for it = 1:numel(nz_r)
        
        if (segFiringPop(nz_r(it), nz_c(it)) == -1)  % new blob
            visitPoints_r = [nz_r(it)];
            visitPoints_c = [nz_c(it)];
            
            while numel(visitPoints_r) ~= 0
                % find neighbours, 8-neighborhood
                r = visitPoints_r(1); visitPoints_r(1) = [];
                c = visitPoints_c(1); visitPoints_c(1) = [];
                r_mod = mod(r - 1, sheet_size) + 1;
                c_mod = mod(c - 1, sheet_size) + 1;
                segFiringPop(r, c) = segn;
                %segFiringPop(r_mod, c_mod) = segn;
                
                neighb_r = [r-1, r-1, r-1, r,   r,   r+1, r+1, r+1];
                neighb_c = [c-1, c,   c+1, c-1, c+1, c-1, c,   c+1];
                
                % wrap around square boundaries (we start indexing at 1)
%                 neighb_r(find(neighb_r < 1)) = sheet_size;
%                 neighb_r(find(neighb_r > sheet_size)) = 1;
%                 neighb_c(find(neighb_c < 1)) = sheet_size;
%                 neighb_c(find(neighb_c > sheet_size)) = 1;
                
                for n_i = 1:8
                    r_mod = mod(neighb_r(n_i) - 1, sheet_size) + 1;
                    c_mod = mod(neighb_c(n_i) - 1, sheet_size) + 1;
                    if (thrFiringPop(r_mod, c_mod) == 1 && ...
                        segFiringPop(neighb_r(n_i), neighb_c(n_i)) == -1)
                        
                        visitPoints_r = [visitPoints_r neighb_r(n_i)];
                        visitPoints_c = [visitPoints_c neighb_c(n_i)];
                        segFiringPop(neighb_r(n_i), neighb_c(n_i)) = segn;
                        %segFiringPop(r_mod, c_mod) = segn;
                        %pcolor(segFiringPop);
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
            segCenters_r(it) = mean(r);
            segCenters_c(it) = mean(c);
            
            % Have all this modulo sheet_size, so when the blob has wrapped, it
            % does not get out of boundaries
            segCenters_r(it) = mod(segCenters_r(it) - 1, sheet_size) + 1;
            segCenters_c(it) = mod(segCenters_c(it) - 1, sheet_size) + 1;
        end
    end
    
    % This is ugly, but the tracking thing produces a mirrored copy of the
    % blob and thus two identical output values
    segCenters_r = segCenters_r(1);
    segCenters_c = segCenters_c(1);
    
end
