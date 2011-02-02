% Compute ISI statistics about neurons in the simulation

%figure;

delete_high_isi = true;
isi_threshold = 0.5; % sec

min_stat = 750; % minimal no.of ISI samples

% load

createISI;  % script to compute ISI cell array - side effects

for nID = 1:sheet_size^2
    tmp_isi = isi{nID};
    if delete_high_isi == true
       tmp_isi(find(tmp_isi > isi_threshold)) = [];
    end
    
    if (numel(tmp_isi) >= min_stat)
        isi_mean(nID) = mean(tmp_isi);
        isi_std(nID) = std(tmp_isi);
    else
        isi_mean(nID) = nan;
        isi_std(nID) = 1;
    end
end

C_v = isi_std./isi_mean;
plot(isi_mean, C_v, 'o', 'MarkerSize', 3);
xlabel('Mean ISI (s)');
ylabel('Coefficient of variation - C_v');