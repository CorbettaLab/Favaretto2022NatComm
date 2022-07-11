%% ===================================================================== %%
%            SUBCORTICAL PCA - CORTICAL/SUBCORTICAL INTERACTION           %
%  ===================================================================== %%

%{
inputs:
- LEs_bold: Leading Eigenvectors
- NET_reduction: output od "dimensionality_reduction.m"
- DFS: output of "main_DFSs.m"
- n: number of ROIs
- n_cort: number of cortical ROIs
- n_net: number of networks
- Time_ctr_pat_bold: output of "eval_LEs_DFC.m"
- n_Subjects: number of subjects in each group

outputs:
- PC_Sub_res: (Struct) with the following fields:
    - SCs_loadings
    - SCs_scores 
    - ExpVar_sub  
    - labels
    - H0_independent: probability of independent shifts of cortical and
    subcortical conenctivity patterns
    - obs_cor_sub_cond_sbj: conditioned probability of subcortical shift,
    given a shift in cortical connectivity pattern
%}

function PC_Sub_res = subcortical_pca(LEs_bold, NET_reduction, DFS, n, n_cort, n_net, Time_ctr_pat_bold, n_Subjects)


% Correlation-based normalization
for ii = 2 : size(LEs_bold, 2)
    if corr(LEs_bold(:, ii-1), LEs_bold(:, ii)) < 0
        LEs_bold(:, ii) = - LEs_bold(:, ii);
    end
end
LEI_subcort = LEs_bold(72 : end, :);


[PC_sub, Scores_sub, ~, ~, ExpVar_sub] = pca(LEI_subcort');
Scores_sub = zscore(Scores_sub);

labels = {'L-CER', 'L-THAL', 'L-CAU', 'L-PUT', 'L-PAL', 'BSTEM', 'L-HIPP', ...
    'L-AMY', 'L-ACC', 'L-DIE', 'R-CER', 'R-THAL', 'R-CAU', 'R-PUT', ...
    'R-PAL', 'R-HIPP', 'R-AMY', 'R-ACC', 'R-DIE'};


PC_Sub_res.SCs_loadings = PC_sub;
PC_Sub_res.SCs_scores = Scores_sub;
PC_Sub_res.ExpVar_sub = ExpVar_sub; 
PC_Sub_res.labels = labels;


%% Project DFS onto PC1 and PC2 

K = size(DFS, 2);

mask        = triu(ones(n), 1);
imask       = find(mask);
[ir, ic]    = ind2sub([n n], imask);

DFS_mat = nan(n, n, K);
DFS_pj  = nan(n_cort, 2, K);    % [cort | PC subcort | K]

DFS = DFS ./ max(abs(DFS), [], 2);
for ii = 1 : length(imask)
    DFS_mat(ir(ii), ic(ii), :) = DFS(:, ii)';
    DFS_mat(ic(ii), ir(ii), :) = DFS(:, ii)'; 
end

for ik = 1 : K
    tmp = squeeze(DFS_mat(1 : n_cort, (n_cort+1) : end, ik));
    DFS_pj(:, :, ik) = (tmp - mean(tmp)) * PC_sub(:, 1 : 2);
end

% Average connectivity within subcortical PCs and cortical Nets
DFS_pj_net  = nan(2, n_net, K);

for in = 1 : n_net 
    idx = find(NET_reduction.net_reduced == in);
    
    DFS_pj_net(:, in, :) = squeeze(nanmean(DFS_pj(idx, :, :), 1));
end


%% Time-course plot

Net_LEI_ave = nan(size(LEs_bold, 2), n_net);

for in = 1 : n_net
    idx                 = find(NET_reduction.net_reduced == in);
    Net_LEI_ave(:, in)  = nanmean(LEs_bold(idx, :), 1)';
end

Net_LEI_ave = zscore(Net_LEI_ave);


%% Time measure reorganization

% Variability 
diff_sub = diff(Scores_sub(:, [1 2]), [], 1);
diff_cor = diff(Net_LEI_ave, [], 1);

thr_diff = 0.2937; % to be derived with elbow method
diff_sub_bin = zeros(size(diff_sub));
diff_sub_bin(abs(diff_sub) > thr_diff) = 1;
diff_cor_bin = zeros(size(diff_cor));
diff_cor_bin(abs(diff_cor) > thr_diff) = 1;
diff_sub_bin_all = double(diff_sub_bin(:, 1) | diff_sub_bin(:, 2));

poiss_prob = @(x, l) l^x * exp(-l) / factorial(x); 

p1_sub      = length(find(diff_sub_bin_all == 1)) ./ length(diff_sub_bin_all); %poiss_prob(1, lambda_sub);
p1_cor      = nan(n_net, 1);
p1_cor_sub  = nan(n_net, 1); % P(cor = 1, sub = 1)
for ii = 1 : n_net
    p1_cor(ii)      = length(find(diff_cor_bin(:, ii) == 1)) ./ ...
        length(diff_cor_bin(:, ii)); % poiss_prob(1, lambda_cor(ii));
    p1_cor_sub(ii)  = p1_cor(ii) * p1_sub;
end

% Since independent: P(sub | cor) = P(sub & cor) / P(cor) = P(sub) * P(cor) / P(cor) = P(sub) 

%% Evaluation of differences subject-wise

Scores_sub_sbj  = split_LEigs_subj(Scores_sub(:, [1 2])', Time_ctr_pat_bold, n_Subjects);
Net_LEI_ave_sbj = split_LEigs_subj(Net_LEI_ave', Time_ctr_pat_bold, n_Subjects); 

Ns              = sum(n_Subjects);
lambda_sub_sbj  = nan(Ns, 1);
lambda_cor_sbj  = nan(Ns, n_net);

p1_sub_sbj      = nan(Ns, 1);
p1_cor_sbj      = nan(Ns, n_net);
p1_cor_sub_sbj  = nan(Ns, n_net);

obs_cor_sub_cond_sbj  = nan(Ns, n_net);   % observed conditioned probability

for is = 1 : Ns
    tmp_sub = squeeze(Scores_sub_sbj(:, :, is));
    tmp_cor = squeeze(Net_LEI_ave_sbj(:, :, is));
    
    tmp_diff_sub = diff(tmp_sub, [], 2);
    tmp_diff_cor = diff(tmp_cor, [], 2);
    
    tmp_diff_sub_bin = zeros(size(tmp_diff_sub));
    tmp_diff_sub_bin(abs(tmp_diff_sub) > thr_diff) = 1;
    tmp_diff_cor_bin = zeros(size(tmp_diff_cor));
    tmp_diff_cor_bin(abs(tmp_diff_cor) > thr_diff) = 1;
    tmp_diff_sub_bin_all = double(tmp_diff_sub_bin(1, :) | tmp_diff_sub_bin(2, :));
    
    lambda_sub_sbj(is)    = poissfit(tmp_diff_sub_bin_all');
    lambda_cor_sbj(is, :) = length(find(tmp_diff_sub_bin_all == 1)) ./ ...
        length(tmp_diff_sub_bin_all); %poissfit(tmp_diff_cor_bin');
    
    p1_sub_sbj(is) = poiss_prob(1, lambda_sub_sbj(is));
    for ii = 1 : n_net
        p1_cor_sbj(is, ii)      = length(find(tmp_diff_cor_bin(ii, :) == 1)) ./ ...
            length(tmp_diff_cor_bin(ii, :)); %poiss_prob(1, lambda_cor_sbj(is, ii));
        p1_cor_sub_sbj(is, ii)  = p1_cor_sbj(is, ii) * p1_sub_sbj(is);
        
        obs_cor_sub_cond_sbj(is, ii) = length(find(tmp_diff_sub_bin_all == 1 & ...
            tmp_diff_cor_bin(ii, :) == 1)) ./ length(find(tmp_diff_cor_bin(ii, :) == 1));
    end
end

%% Wilcoxon rank sum test

for in = 1 : n_net
    [p, h, st] = ranksum(obs_cor_sub_cond_sbj(:, in), p1_sub_sbj);
    disp(['Net ' , lab{in}, ':'])
    disp(['z = ', num2str(st.zval), ' (p = ', num2str(p), '; pbonf = ', ...
        num2str(p * n_net), ')'])
end

PC_Sub_res.H0_independent = p1_cor_sub_sbj;
PC_Sub_res.obs_cor_sub_cond_sbj = obs_cor_sub_cond_sbj;





