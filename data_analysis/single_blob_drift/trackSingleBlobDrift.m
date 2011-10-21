function [blobPos_r, blobPos_c] = trackSingleBlobDrift(startTime, endTime, spikeHist, dt_track, delta_t)

    sheet_size = sqrt(size(spikeHist, 1));
    minEdgeDist =  10;

    blobPos_r = [0];
    blobPos_c = [0];

    t_start = startTime;
    t_end = endTime-delta_t/2;

    firingPop = getFiringPop(spikeHist, t_start, dt_track, delta_t);
    %firingPop = reshape(firingPop, sheet_size, sheet_size)';
    [currBlobPos_r currBlobPos_c] = trackSingleBlob(firingPop);

    shift_r = currBlobPos_r - blobPos_r(end);
    shift_c = currBlobPos_c - blobPos_c(end);

    updateId = 0;
    nUpdate = 100;
    for t = t_start:dt_track:t_end
        firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
        %firingPop = reshape(firingPop, sheet_size, sheet_size)';
        [new_r, new_c] = trackBlobs(firingPop);

        % Check if the blob has been wrapped around boundary
        if (abs(new_r - blobPos_r(end)) > sheet_size/2 || ...
            abs(new_c - blobPos_c(end)) > sheet_size/2)
        
            disp('Blob wrap-around');
            shift_r = currBlobPos_r - blobPos_r(end);
            shift_c = currBlobPos_c - blobPos_c(end);
        else
            currBlobPos_r = new_r;
            currBlobPos_c = new_c;
        end
        
        blobPos_r = [blobPos_r (r - shift_r)];
        blobPos_c = [blobPos_c (c - shift_c)];



        if (mod(updateId, nUpdate) == 0)
            t
        end
        updateId = updateId + 1;
    end
end