%% ===================================================================== %%
%                          Dynamic measures PCA                           %
%  ===================================================================== %%
% 
% inputs:
% - Dyn_all: [n_sub, K*(K+1)] matrix of all dynamical measures (output of "DFSs_dynamics.m")
% - K: number of DFSs
% - idx_subj_cond_st: [1, sum(n_Subjects)]: vector to associate each subject to a group: 
%      - 1: CTR, 
%      - 2.1: PAT(2w-severe), 2.2: PAT(2w-mild)
%      - 3.1: PAT(3m-severe), 3.2: PAT(3m-mild)
%      - 4.1: PAT(1y-severe), 4.2: PAT(1y-mild)
% 
% outputs:
% - Dyn_PCs: (Struct)
%     - Dyn_PCs.PCs_dyn: loading of Dyn-PCs;
%     - Dyn_PCs.Scores_dyn = scores of Dyn-PCs (CTR);
%     - Dyn_PCs.ExpVar_dyn = explained variance;
%     - Dyn_PCs.Scores_xx = scores of Dyn-PCs (PAT)
%     - Dyn_PCs.nPC: number of selected Dyn-PCs
%     - Dyn_PCs.Corr_mat_sort = correlation matrix among dynamical measures
%     - Dyn_PCs.dyn_lab_sort = labels for the correlation matrix
%  ===================================================================== %%


function Dyn_PCs = dynamic_measures_PCA(Dyn_all, K, idx_subj_cond_st)

%% Create Dynamical labels

dyn_lab = cell(size(Dyn_all, 2), 1);

mask_K      = ones(K, K) - eye(K, K);
mask_K_ind  = find(mask_K);
TR_lab = cell(K, K);

for ik = 1 : K
    for iik = 1 : K
        if ik ~= iik
            TR_lab{ik, iik} = ['DFS', num2str(ik), '$>$', num2str(iik)];
        end
    end
end
TR_lab = TR_lab(mask_K_ind);

for ik = 1 : K
    dyn_lab{ik}     = ['$f_', num2str(ik), '$'];
    dyn_lab{ik+K}   = ['$\ell_', num2str(ik), '$'];
end

dyn_lab(2*K+1 : end) = TR_lab;

%% PCA on dynamical data (only controls)

Dyn_ctr = Dyn_all(idx_subj_cond_st == 1, :);
Dyn_all_zCtr = (Dyn_all - nanmean(Dyn_ctr, 1)) ./ nanstd(Dyn_ctr, [], 1);

[PCs_dyn, Scores_dyn, Eig2_dyn, ~, ExpVar_dyn, mu_dyn] = ...
    pca(Dyn_all_zCtr(idx_subj_cond_st == 1, :));

Ns = length(find(idx_subj_cond_st > 1 & idx_subj_cond_st < 3)); % number of patients
Scores_pat = nan(3*Ns, size(Scores_dyn, 2));
for is = 1 : (3*Ns)
    data = Dyn_all_zCtr(is + length(find(idx_subj_cond_st == 1)), :);
    Scores_pat(is, :) = data * PCs_dyn;
end


% PCs scaling
for ip = 1 : size(PCs_dyn, 2)
    max_val             = max(abs(PCs_dyn(:, ip)));
    PCs_dyn(:, ip)      = PCs_dyn(:, ip) ./ max_val;
    Scores_dyn(:, ip)   = Scores_dyn(:, ip) .* max_val;
    Scores_pat(:, ip)   = Scores_pat (:, ip) .* max_val;
end

% Reproject patients' dynamical measures onto the PCs space
Scores_pat = (Scores_pat - nanmean(Scores_dyn)) ./ nanstd(Scores_dyn);

Scores_2w = Scores_pat(1 : Ns, :);
Scores_3m = Scores_pat(Ns+1 : 2*Ns, :);
Scores_1y = Scores_pat(2*Ns+1 : 3*Ns, :);

%% Save data

Dyn_PCs.PCs_dyn = PCs_dyn;
Dyn_PCs.Scores_dyn = Scores_dyn;
Dyn_PCs.ExpVar_dyn = ExpVar_dyn;
Dyn_PCs.Scores_2w = Scores_2w;
Dyn_PCs.Scores_3m = Scores_3m;
Dyn_PCs.Scores_1y = Scores_1y;

    

idx_eig = find(Eig2_dyn > 1);
idx_var = find(ExpVar_dyn > 10);

nPC = min([length(idx_eig), length(idx_var)]);
Dyn_PCs.nPC = nPC;

%% Correlation of Dynamical measures (order w.r.t. PCs)

[~, Cl] = max(abs(PCs_dyn(:, 1 : nPC)), [], 2);
idx_cc = [];
for ii = 1 : nPC
    idx_pc  = find(Cl == ii);
    PCs_tmp = PCs_dyn(idx_pc, ii);
    [~, isort] = sort(PCs_tmp, 'descend');
    idx_cc = [idx_cc; idx_pc(isort)];
end

CC = corr(Dyn_all_zCtr(idx_subj_cond_st == 1, :));
CC = CC(idx_cc, idx_cc);

Dyn_PCs.Corr_mat_sort = CC;
Dyn_PCs.dyn_lab_sort = dyn_lab(idx_cc);


end




    
    
