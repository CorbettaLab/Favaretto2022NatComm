%% Posthoc permutation test for Transition probabilities
% 
% inputs:
% - Trans_all: [K(K-1), n_sub]
% - K: number of DFSs
% - idx_subj_cond_st: (vector [n_sub, 1]) integer indicating the group 
%     of each subject: 1: CTR, 
%                      2.1: PAT(2w-severe), 2.2: PAT(2w-mild)
%                      3.1: PAT(3m-severe), 3.2: PAT(3m-mild)
%                      4.1: PAT(1y-severe), 4.2: PAT(1y-mild)

function [pval_Trans, pval_Trans_fdr, pval_Trans_posthoc] = ...
    posthoc_Trans_ctr_pat(Trans_all, K, idx_subj_cond_st)

addpath('../utility/')

pval_Trans = nan(K, K);

for ik = 1 : K
    for iik = setdiff(1 : K, ik)
        data = squeeze(Trans_all(ik, iik, idx_subj_cond_st < 3));
        pval_Trans(ik, iik) = ...
            kruskalwallis(data, idx_subj_cond_st(idx_subj_cond_st < 3), 'off');
    end
end


% FDR Correction
idx_T = find(~isnan(pval_Trans));
[r_tmp, c_tmp] = ind2sub([K K], idx_T);
pval_fdr = mafdr(pval_Trans(idx_T), 'BHFDR', true);

pval_Trans_fdr = nan(K, K);

for ii = 1 : length(idx_T)
    pval_Trans_fdr(r_tmp(ii), c_tmp(ii)) = pval_fdr(ii);
end

pval_Trans_posthoc = cell(K, K, 3); 

for ik = 1 : K
    for iik = setdiff(1 : K, ik)
        
        % if significant: post-hoc 
        if pval_Trans_fdr(ik, iik) < .05
            
            % CTRs vs PAT2w
            cond_tmp = unique(idx_subj_cond_st(idx_subj_cond_st < 3));
            cond_tmp(isnan(cond_tmp)) = [];
            pval_tmp = nan(length(cond_tmp));
            m_tmp = nan(length(cond_tmp), 1);
            s_tmp = nan(length(cond_tmp), 1);
            for ic = 1 : length(cond_tmp)
                a = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(ic)));
                m_tmp(ic) = nanmean(a);
                s_tmp(ic) = nanstd(a) ./ sqrt(length(a));
                for iic = (ic+1) : length(cond_tmp)
                    b = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(iic)));
                    stats = permutation_htest2_np([a; b]', ...
                        [ones(1, numel(a)) 2*ones(1, numel(b))], 2000, 0.05, 'ttest');

                    pval_tmp(ic, iic) = min(stats.pvals);
                    pval_tmp(iic, ic) = pval_tmp(ic, iic);
                end
            end
            pval_tmp(pval_tmp > .05) = nan;
            pval_Trans_posthoc{ik, iik, 1} = pval_tmp;
            
            
            % Check what happens at 3 months
            cond_tmp = unique(idx_subj_cond_st(idx_subj_cond_st < 2 | ...
                idx_subj_cond_st == 3.1 | idx_subj_cond_st == 3.2));
            cond_tmp(isnan(cond_tmp)) = [];
            pval_tmp = nan(length(cond_tmp));
            m_tmp = nan(length(cond_tmp), 1);
            s_tmp = nan(length(cond_tmp), 1);
            for ic = 1 : length(cond_tmp)
                a = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(ic)));
                m_tmp(ic) = nanmean(a);
                s_tmp(ic) = nanstd(a) ./ sqrt(length(a));
                for iic = (ic+1) : length(cond_tmp)
                    b = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(iic)));
                    stats = permutation_htest2_np([a; b]', ...
                        [ones(1, numel(a)) 2*ones(1, numel(b))], 2000, 0.05, 'ttest');

                    pval_tmp(ic, iic) = min(stats.pvals);
                    pval_tmp(iic, ic) = pval_tmp(ic, iic);
                end
            end
            pval_tmp(pval_tmp > .05) = nan;
            pval_Trans_posthoc{ik, iik, 2} = pval_tmp;
            
            
            % Check what happens at 1 year
            cond_tmp = unique(idx_subj_cond_st(idx_subj_cond_st < 2 | ...
                idx_subj_cond_st > 3));
            cond_tmp(isnan(cond_tmp)) = [];
            pval_tmp = nan(length(cond_tmp));
            m_tmp = nan(length(cond_tmp), 1);
            s_tmp = nan(length(cond_tmp), 1);
            for ic = 1 : length(cond_tmp)
                a = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(ic)));
                m_tmp(ic) = nanmean(a);
                s_tmp(ic) = nanstd(a) ./ sqrt(length(a));
                for iic = (ic+1) : length(cond_tmp)
                    b = squeeze(Trans_all(ik, iik, idx_subj_cond_st == cond_tmp(iic)));
                    stats = permutation_htest2_np([a; b]', ...
                        [ones(1, numel(a)) 2*ones(1, numel(b))], 2000, 0.05, 'ttest');

                    pval_tmp(ic, iic) = min(stats.pvals);
                    pval_tmp(iic, ic) = pval_tmp(ic, iic);
                end
            end
            pval_tmp(pval_tmp > .05) = nan;
            pval_Trans_posthoc{ik, iik, 3} = pval_tmp;

        end
    end
end

rmpath('../utility/')
