function [Time_ctr_pat_bold, LEs_bold, Eigs_bold, DFC_bold] = ...
    bold_DFC_LEs_extraction(n_Subjects, n_Cond, Time_ctr_pat, n, Tmin, win_len, i_mask, ...
    l_mask, BOLD_ctr_pat)

LEs_bold = zeros(n, (Tmin - win_len) * sum(n_Subjects));
Eigs_bold  = zeros(n, (Tmin - win_len) * sum(n_Subjects));
ik = 1;

DFC_bold     = zeros(l_mask, (Tmin - win_len) * sum(n_Subjects));
iik = 1;


for cond = 1 : (n_Cond + 1)
    
    disp(['CONDITION ', num2str(cond)])
    for is = 1 : n_Subjects(cond)
        
        disp(['Subject ', num2str(is)])
        idx = find((Time_ctr_pat(1, :) == cond) & (Time_ctr_pat(2, :) == is));
        
        n_win   = length(idx) - win_len;   % # sliding windows
        DFC     = zeros(n_win, l_mask);
        BOLD_is = BOLD_ctr_pat(idx, :);
        
        for it = 1 : n_win
            data          = BOLD_is(it : it + (win_len-1), :)';
            
            % istantaneous FC
            FC            = atanh(corr(data', 'rows', 'pairwise'));
            FC(isnan(FC)) = 0;
            FC(isinf(FC)) = 0;
            
            DFC(it, :) = FC(i_mask);
            
            [V1, E1] = eigs(FC, n);
            LEs_bold(:, ik) = V1(:, 1);
            Eigs_bold(:, ik)     = diag(abs(E1));
            
            ik = ik + 1;
            
            DFC_bold(:, iik) = DFC(it, :);
            
            iik = iik + 1;
        end
    end
    
    
end

l_bold              = size(Time_ctr_pat, 2) - sum(win_len .* n_Subjects);
Time_ctr_pat_bold   = zeros(3, l_bold);

ik = 0;
for cond = 1 : (n_Cond + 1)
    for is = 1 : n_Subjects(cond)
        idx_T = find((Time_ctr_pat(1, :) == cond) & (Time_ctr_pat(2, :) == is));
        
        Time_ctr_pat_bold(:, ik+1 : ik+(Tmin-win_len)) = ...
            Time_ctr_pat(:, idx_T(1 : Tmin-win_len));
        
        ik = ik + (Tmin - win_len);
    end
end
