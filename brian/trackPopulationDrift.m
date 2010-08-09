function [blobPos_r, blobPos_c] = trackPopulationDrift(startTime, endTime, spikeHist, dt_track, delta_t)

    sheet_size = sqrt(size(spikeHist, 1));
    minEdgeDist =  10;

    blobPos_r = [0];
    blobPos_c = [0];

    t_start = startTime;
    t_end = endTime-delta_t/2;

    firingPop = getFiringPop(spikeHist, t_start, dt_track, delta_t);
    [currBlobPos_r currBlobPos_c] = trackBlobs(firingPop);

    % Find a blob nearest to the center of sheet
    c_r = sheet_size/2;
    c_c = sheet_size/2;
    [cMin cMin_i] = min((currBlobPos_r - c_r).^2 + (currBlobPos_c - c_c).^2);

    currBlobPos_r = currBlobPos_r(cMin_i);
    currBlobPos_c = currBlobPos_c(cMin_i);

    shift_r = currBlobPos_r - blobPos_r(end);
    shift_c = currBlobPos_c - blobPos_c(end);

    updateId = 0;
    nUpdate = 100;
    for t = t_start:dt_track:t_end
        firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
        [r, c] = trackBlobs(firingPop);

        % Find blob nearest to the last position and record its position
        dist = (r-currBlobPos_r).^2 + (c-currBlobPos_c).^2;
        [minDist min_i] = min(dist);

        new_r = r(min_i);
        new_c = c(min_i);

        blobPos_r = [blobPos_r (new_r - shift_r)];
        blobPos_c = [blobPos_c (new_c - shift_c)];

        % Check if new position of the blob is too near the edge, if so
        % choose new blob near the center and reset shift
        if (sheet_size - new_r < minEdgeDist || new_r < minEdgeDist || ...
            sheet_size - new_c < minEdgeDist || new_c < minEdgeDist)

            disp('change of blob');
            [cMin cMin_i] = min((r - c_r).^2 + (c - c_c).^2);
            currBlobPos_r = r(cMin_i);
            currBlobPos_c = c(cMin_i);

            shift_r = currBlobPos_r - blobPos_r(end);
            shift_c = currBlobPos_c - blobPos_c(end);
        else
            currBlobPos_r = new_r;
            currBlobPos_c = new_c;
        end

        if (mod(updateId, nUpdate) == 0)
            t
        end
        updateId = updateId + 1;
    end
end