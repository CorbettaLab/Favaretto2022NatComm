%% Posthoc permutation test for fraction times or dwell time
% 
% inputs:
% - data_ctr [K, n_CTR]
% - data_pat [K, n_PAT(2w, 3m or 1y)]
% - ST [n_PAT, 1] 1: severe, 0: mild


function [pval_ctr_pat, pval_ctr_pat_fdr] = posthoc_ctr_pat_perm(data_ctr, data_pat, ST)

addpath('../utility/')

K = size(data_pat, 1);
st_severe = find(ST == 1);
st_mild = find(ST == 0);

pval_ctr_pat = nan(3*K, 3*K);

for ik = 1 : K
    a = data_ctr(ik, :);
    b = data_pat(ik, st_severe);
    c = data_pat(ik, st_mild);
    
    % a vs b
    stats = permutation_htest2_np([a, b], ...
        [ones(1, numel(a)) 2*ones(1, numel(b))], 2000, 0.05, 'ttest');
    pval_ctr_pat(3*(ik-1)+1, 3*(ik-1)+2) = min(stats.pvals);
    
    % a vs c
    stats = permutation_htest2_np([a, c], ...
        [ones(1, numel(a)) 2*ones(1, numel(c))], 2000, 0.05, 'ttest');
    pval_ctr_pat(3*(ik-1)+1, 3*(ik-1)+3) = min(stats.pvals);
    
    % b vs c
    stats = permutation_htest2_np([b, c], ...
        [ones(1, numel(b)) 2*ones(1, numel(c))], 2000, 0.05, 'ttest');
    pval_ctr_pat(3*(ik-1)+2, 3*(ik-1)+3) = min(stats.pvals);
end

idx         = find(~isnan(pval_ctr_pat));
[r, c]      = ind2sub(size(pval_ctr_pat), idx);
p_tmp       = pval_ctr_pat(idx);
pfdr_tmp    = mafdr(p_tmp, 'BHFDR', true);

pval_ctr_pat_fdr = nan(size(pval_ctr_pat));
for ii = 1  : length(idx)
    pval_ctr_pat_fdr(r(ii), c(ii)) = pfdr_tmp(ii);
    pval_ctr_pat_fdr(c(ii), r(ii)) = pfdr_tmp(ii);
end

pval_ctr_pat_fdr(pval_ctr_pat_fdr > 0.05) = nan;

mask_pos = nan(3*K, 3*K);
mask_neg = nan(3*K, 3*K);
for ik = 1 : K
    mask_pos((ik-1)*3+1, (ik-1)*3+2) = 1;
    mask_neg((ik-1)*3+1, (ik-1)*3+3) = 1;
end

mask_pos = pval_ctr_pat_fdr .* mask_pos;
mask_neg = pval_ctr_pat_fdr .* mask_neg;

mask_pos = mask_pos(1 : 3 : end, 2 : 3 : end);
mask_neg = mask_neg(1 : 3 : end, 3 : 3 : end);

pval_ctr_pat_fdr = [nanmin([diag(mask_pos), diag(mask_neg)], [], 2)';
    diag(mask_pos)'; diag(mask_neg)'];

pval_ctr_pat_fdr(pval_ctr_pat_fdr == 0) = nan;

rmpath('../utility/')