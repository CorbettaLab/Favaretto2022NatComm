function LEigs_time = split_LEigs_subj(LEigs, Time_ctr_pat, n_Subjects)

Ns  	= sum(n_Subjects);  % # subjects (total)
[n, T] 	= size(LEigs);      % [# nodes, # time points]

n_Cond = length(n_Subjects);  % # conditions

nT = T / Ns;    % length of each time series

LEigs_time = zeros(n, nT, Ns);

    
iis = 1;
for cond = 1 : n_Cond
    for is = 1 : n_Subjects(cond)
        indT = (Time_ctr_pat(1, :) == cond) & ...
            (Time_ctr_pat(2, :) == is);
        
        LEigs_time(:, :, iis) = LEigs(:, indT);
        
        iis = iis + 1;
    end
end



            