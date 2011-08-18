function fig_Cv_stat(res, min_stat, nBins, fontSize)
    % ISI statistics
    % Compute ISI statistics about neurons in the simulation


    isi_vec{1} =  createISI(res.spikeCell_e);
    isi_vec{2} =  createISI(res.spikeCell_i);
    tit{1} = 'Principal cells';
    tit{2} = 'Interneurons';

    for it = 1:2
        isi = isi_vec{it};
    
        N = size(isi, 2);

        for nID = 1:N
            tmp_isi = isi{nID};

            if (numel(tmp_isi) >= min_stat)
                isi_mean(nID) = mean(tmp_isi);
                isi_std(nID) = std(tmp_isi);
            else
                isi_mean(nID) = nan;
                isi_std(nID) = 1;
            end
        end

        subplot(1, 2, it, 'FontSize', fontSize);
        C_v = isi_std./isi_mean;
        [hist_N, hist_X] = hist(C_v, nBins);
        bar(hist_X, hist_N/N);
        xlabel('C_v');
        ylabel('No. of cells');
        title(tit{it});
        %xlim([0 2]);
        
        clear isi_mean isi_std isi
    end
end