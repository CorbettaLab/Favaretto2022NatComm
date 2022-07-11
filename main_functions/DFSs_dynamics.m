%% ===================================================================== %%
%                  DYNAMICS of DYNAMICAL FUNCTIONAL STATES                %
%  ===================================================================== %%
% 
% inputs:
% - Time_ctr_pat_bold: [3, nW] (output of "eval_LEs_DFS.m")
% - n_Subjects: [1, 4] number of subjects in each group
% - K: number of DFSs
% - DFS_sw: [nW, 1] (output of "main_DFS.m"))
% - Static_info; (Struct) (output of "static_severity_evaluation.m")
% - idx_subj_cond [1, sum(n_Subjects)] (output of "eval_LEs_DFC.m")
% 
% outputs:
% - Dyn_all: horizontal concatenation of:
%     - fraction_time [sum(n_Subjects), K]
%     - dwell_time [sum(n_Subjects), K]
%     - Trans_prob (vectorized) [sum(n_Subjects), K*(K-1)]
% - idx_subj_cond_st [1, sum(n_Subjects)]: as idx_subj_cond with additional
% distinction between severe and mild patients:
%     1: CTR, 
%     2.1: PAT(2w-severe), 2.2: PAT(2w-mild)
%     3.1: PAT(3m-severe), 3.2: PAT(3m-mild)
%     4.1: PAT(1y-severe), 4.2: PAT(1y-mild)
%  ===================================================================== %%

function [Dyn_all, idx_subj_cond_st] = DFSs_dynamics(Time_ctr_pat_bold, n_Subjects, K, DFS_sw, ...
    Static_info, idx_subj_cond)

%% Add functions to path
addpath(genpath('utility/'))
addpath('external/BCT')

Tctr_bold = find(Time_ctr_pat_bold(1, :) == 1);
T2w_bold  = find(Time_ctr_pat_bold(1, :) == 2);
T3m_bold  = find(Time_ctr_pat_bold(1, :) == 3);
T1y_bold  = find(Time_ctr_pat_bold(1, :) == 4);

%% Patients Split (severe - mild)
idx_subj_cond_st = idx_subj_cond;
st_severe = find(Static_info.ST > 0);
st_mild = find(Static_info.ST < 0);

idx1 = find(idx_subj_cond == 2);
idx_subj_cond_st(idx1(st_severe)) = 2.1;
idx_subj_cond_st(idx1(st_mild)) = 2.2;
idx1 = find(idx_subj_cond == 3);
idx_subj_cond_st(idx1(st_severe)) = 3.1;
idx_subj_cond_st(idx1(st_mild)) = 3.2;
idx1 = find(idx_subj_cond == 4);
idx_subj_cond_st(idx1(st_severe)) = 4.1;
idx_subj_cond_st(idx1(st_mild)) = 4.2;

%% Dynamics of DFSs

% Split DFS_sw (= DFS sequence for each subject)
DFS_sw_sbj = squeeze(split_LEigs_subj(DFS_sw', Time_ctr_pat_bold, n_Subjects));

%% Evaluate DFSs' frequency of occurrence

frac_dfs_ctr = crosstab(DFS_sw(Tctr_bold), Time_ctr_pat_bold(2, Tctr_bold));
frac_dfs_ctr = frac_dfs_ctr ./ sum(frac_dfs_ctr, 1);

frac_dfs_pats_2w = crosstab(DFS_sw(T2w_bold), Time_ctr_pat_bold(2, T2w_bold));
frac_dfs_pats_2w = frac_dfs_pats_2w ./ sum(frac_dfs_pats_2w, 1);

frac_dfs_pats_3m = crosstab(DFS_sw(T3m_bold), Time_ctr_pat_bold(2, T3m_bold));
frac_dfs_pats_3m = frac_dfs_pats_3m ./ sum(frac_dfs_pats_3m, 1);

frac_dfs_pats_1y = crosstab(DFS_sw(T1y_bold), Time_ctr_pat_bold(2, T1y_bold));
frac_dfs_pats_1y = frac_dfs_pats_1y ./ sum(frac_dfs_pats_1y, 1);

frac_all = cat(2, frac_dfs_ctr, frac_dfs_pats_2w, frac_dfs_pats_3m, ...
    frac_dfs_pats_1y);


%% States Dwell Time distributions

% Subject-wise Dwell-Time

DT_sbj  = nan(sum(n_Subjects), K);

for is = 1 : sum(n_Subjects)
    DFS_is = DFS_sw_sbj(:, is);
    dwell_tmp    = evaluate_dwelltime(DFS_is, K);
    DT_sbj(is, :)   = cellfun(@(x) mean(x), dwell_tmp);
end
DT_sbj(isnan(DT_sbj)) = 0;


%% Evaluate States Transition

n_CTR       = n_Subjects(1);
n_PAT_all   = n_Subjects(2);

Trans_ctr       = nan(K, K, n_CTR);
Trans_ctr_v2    = nan(K, K, n_CTR);
Trans_pat       = {nan(K, K, n_PAT_all); nan(K, K, n_PAT_all); nan(K, K, n_PAT_all)};
Trans_pat_v2    = Trans_pat;

Jumps_ctr = nan(K, K, n_CTR);
Jumps_pat = {nan(K, K, n_PAT_all); nan(K, K, n_PAT_all); nan(K, K, n_PAT_all)};

% CTR
DFS_tmp     = DFS_sw(Tctr_bold);
Time_tmp    = Time_ctr_pat_bold(:, Tctr_bold);
for is = 1 : n_CTR
    idx                 = find(Time_tmp(2, :) == is);
    [Trans_ctr(:, :, is), idx_trans, Trans_ctr_v2(:, :, is)] = ...
        evaluate_DFS_transition_noDiagonal(DFS_tmp(idx), K);
    Jumps_ctr(:, :, is) = evaluate_DFS_jumps(DFS_tmp(idx), K);
end


% PAT
for cond = 2 : 4
    switch cond
        case 2
            DFS_tmp     = DFS_sw(T2w_bold);
            Time_tmp    = Time_ctr_pat_bold(:, T2w_bold);
        case 3
            DFS_tmp     = DFS_sw(T3m_bold);
            Time_tmp    = Time_ctr_pat_bold(:, T3m_bold);
        case 4
            DFS_tmp     = DFS_sw(T1y_bold);
            Time_tmp    = Time_ctr_pat_bold(:, T1y_bold);
    end
    for is = 1 : n_PAT_all
        idx                         = find(Time_tmp(2, :) == is);
        [Trans_pat{cond-1}(:, :, is), idx_trans, Trans_pat_v2{cond-1}(:, :, is)] = ...
            evaluate_DFS_transition_noDiagonal(DFS_tmp(idx), K);
        Jumps_pat{cond-1}(:, :, is) = evaluate_DFS_jumps(DFS_tmp(idx), K);
    end
end


%% Transition vectors

Trans_all  = cat(3, Trans_ctr, Trans_pat{1}, Trans_pat{2}, Trans_pat{3});

mask_K      = ones(K, K) - eye(K, K);
mask_K_ind  = find(mask_K);

Trans_all_vec = nan(length(mask_K_ind), sum(n_Subjects));

for is = 1 : sum(n_Subjects)
    Ttmp = squeeze(Trans_all(:, :, is));
    Trans_all_vec(:, is) = Ttmp(mask_K_ind);
end


%% Put all dynamical information together

Dyn_all = [frac_all', DT_sbj, Trans_all_vec'];

%% REMOVE TOOLBOXES FROM PATH
rmpath(genpath('utility/'))
rmpath('external/BCT')


