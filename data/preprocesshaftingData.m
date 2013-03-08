% Preprocess hafting et al., 2005 matlab data in order to remove NaN values
% and make the time spacing between the rat positions uniform
clear all;

dir = 'hafting_et_al_2005/';
files = [
    'Hafting_Fig2c_Trial1_preprocessed.mat';
    'Hafting_Fig2c_Trial2_preprocessed.mat';
    'Hafting_Fig2d_Trial1_preprocessed.mat';
    'Hafting_Fig2d_Trial2_preprocessed.mat'];
numFiles = [1];
    
tolerance = 1e-6;
tstep = 0.02;

insStart = true;
doCorrection = true;

for f_ind = numFiles;
    fileName = [dir files(f_ind, :)]
    load(fileName);

    if insStart
        % Insert zero velocity inputs at the beginning of position data
        ins_time = 5; % sec
        pos_x = [pos_x(1); pos_x];
        pos_y = [pos_y(1); pos_y];
        pos_timeStamps = pos_timeStamps + ins_time;
        pos_timeStamps = [0; pos_timeStamps];
    end
    
    badIndices = [];

    for i = 1:numel(pos_timeStamps) - 1
        diff = pos_timeStamps(i+1) - pos_timeStamps(i);
        tol = abs(diff - 0.02);
        if (tol > tolerance)
            badIndices = [badIndices i];
        end
    end
    
    if doCorrection
        new_pos_x = [];
        new_pos_y = [];
        new_pos_timeStamps = [];
        
        pos_i = 1;
        for ind = badIndices
            while pos_i <= ind
                new_pos_x = [new_pos_x pos_x(pos_i)];
                new_pos_y = [new_pos_y pos_y(pos_i)];
                new_pos_timeStamps = [new_pos_timeStamps pos_timeStamps(pos_i)];
                
                pos_i = pos_i + 1;
            end
            
            insert = [];
            last = pos_timeStamps(ind) + tstep;
            while last < pos_timeStamps(ind+1)
                if pos_timeStamps(ind+1) - last > tolerance
                    insert = [insert last];
                end
                last = last + tstep;
            end
            
            new_pos_timeStamps = [new_pos_timeStamps insert];
            new_pos_x = [new_pos_x (zeros(1, numel(insert))+pos_x(ind))];
            new_pos_y = [new_pos_y (zeros(1, numel(insert))+pos_y(ind))];
        end
        
        % fill the last part after the last bad index
        while pos_i <= numel(pos_timeStamps)
            new_pos_x = [new_pos_x pos_x(pos_i)];
            new_pos_y = [new_pos_y pos_y(pos_i)];
            new_pos_timeStamps = [new_pos_timeStamps pos_timeStamps(pos_i)];

            pos_i = pos_i + 1;
        end
        
         pos_x = new_pos_x';
         pos_y = new_pos_y';
         pos_timeStamps = new_pos_timeStamps';
    end
    
    save(fileName, 'pos_timeStamps', 'pos_x', 'pos_y');
end
