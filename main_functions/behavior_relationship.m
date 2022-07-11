%% ===================================================================== %%
%             Relationship between behavior and dynamical measures        %                         
%  ===================================================================== %%
% 
% inputs:
% - folder_behavior: folder that contains the following file:
%     - behavior_factor.mat (including 'beh_pat_zCtr', 'task_list')
%     - NIH_tot.mat (including 'NIH_2w', 'NIH_3m', 'NIH_1y')
% - idx_PAT_all: index of patients used for the dynamical analysis
% - Static_info: (Struct) static FC information (output of "static_severity_evaluation.m")
% - Dyn_PCs: (Struct) results of Dynamical measures PCA (output of
% "dynamical_measures_PCA.m"
% 
% outputs:
% - res_pat2w: (Struct) results of the Likelihood ratio test used to determine 
% how Dynamic and Static information explain behavioral deficits at the sub-acute stage 
% - res_recovery: (Struct) results of the Likelihood ratio test used to determine 
% how Dynamic and Static information explain behavioral recovery after 1 year 
%  ===================================================================== %%


function [res_pat2w, res_recovery] =  behavior_relationship(folder_behavior, ...
    idx_PAT_all, Static_info, Dyn_PCs)


load([folder_behavior, '/behavior_factor.mat'], 'beh_pat_zCtr', 'task_list')
load([folder_behavior, '/NIH_tot.mat'], 'NIH_2w', 'NIH_3m', 'NIH_1y')

for it = 1 : length(task_list)
    name = char(task_list{it});
    idx = find(name == '_');
    if ~isempty(idx)
        name(idx) = '-';
        task_list{it} = name;
    end
end

[~, idx] = intersect(NIH_2w.id, idx_PAT_all);
nih_2w = NIH_2w.nih_total(idx);
[~, idx] = intersect(NIH_3m.id, idx_PAT_all);
nih_3m = NIH_3m.nih_total(idx);
[~, idx] = intersect(NIH_1y.id, idx_PAT_all);
nih_1y = NIH_1y.nih_total(idx);


% Change VF and AttValDis
beh_pat_zCtr(:, :, [4 6]) = - abs(beh_pat_zCtr(:, :, [4 6]));

beh_pat_zCtr = beh_pat_zCtr(:, idx_PAT_all, :);
beh_2w_zCtr  = squeeze(beh_pat_zCtr(1, :, :));
beh_3m_zCtr  = squeeze(beh_pat_zCtr(2, :, :));
beh_1y_zCtr  = squeeze(beh_pat_zCtr(3, :, :));

beh_all_pat_zCtr    = cat(1, beh_2w_zCtr, beh_3m_zCtr, beh_1y_zCtr);

% Evaluate Behavioral PCs zCtr (no Visual)
% add NIH
beh_all_pat_zCtr = [beh_all_pat_zCtr, [nih_2w; nih_3m; nih_1y]];
zBeh_zCtr = (beh_all_pat_zCtr - nanmean(beh_all_pat_zCtr)) ./ ...
    nanstd(beh_all_pat_zCtr);



%% Test if Dyn(scores)+FC is better than FC

nPC = Dyn_PCs.nPC;


% Likelihood ratio test of model specification
chi2_ll_pc  = nan(size(zBeh_zCtr, 2), 2);
pvals_ll_pc = nan(size(zBeh_zCtr, 2), 2);
R2_pc       = nan(size(zBeh_zCtr, 2), 3); % [reduced (static) | reduced (dyn) | unreduced]

Ns = length(nih_2w);

for it = 1 : length(pvals_ll_pc)
    
    Y = zBeh_zCtr(1 : Ns, it);
        
    % Restricted models (static)
    [b, dev, stats] = glmfit(Static_info.ST, Y);
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    rL1             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_r1 = glmval(b, Static_info.ST, 'identity');
    R2_pc(it, 1) = 1 - nansum((Y - yhat_r1) .^2) / nansum((Y - nanmean(Y)).^2);
    
    % Restricted models (dynamic)
    [b, dev, stats] = glmfit(Dyn_PCs.Scores_2w(:, 1 : nPC), Y);
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    rL2             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_r2 = glmval(b, Dyn_PCs.Scores_2w(:, 1 : nPC), 'identity');
    R2_pc(it, 2) = 1 - nansum((Y - yhat_r2) .^2) / nansum((Y - nanmean(Y)).^2);
    
    % Unrestricted model
    [b, dev, stats] = glmfit([Dyn_PCs.Scores_2w(:, 1 : nPC), ...
        Static_info.ST], zBeh_zCtr(1 : Ns, it));
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    urL             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_ur = glmval(b, [Dyn_PCs.Scores_2w(:, 1 : nPC), ...
        Static_info.ST], 'identity');
    R2_pc(it, 3) = 1 - nansum((Y - yhat_ur) .^2) / nansum((Y - nanmean(Y)).^2);
    
    [~, pvals_ll_pc(it, 1), chi2_ll_pc(it, 1)] = lratiotest(urL, rL1, nPC);
    [~, pvals_ll_pc(it, 2), chi2_ll_pc(it, 2)] = lratiotest(urL, rL2, 1);

end

res_pat2w.chi2_ll_pc = chi2_ll_pc;
res_pat2w.pvals_ll_pc = pvals_ll_pc; 
res_pat2w.R2_pc = R2_pc;
res_pat2w.task_list = task_list;

%% Plot Dyn measures (scores) vs Recovery Behavior correlation

% Likelihood ratio test of model specification
chi2_ll_1y_pc = nan(size(zBeh_zCtr, 2), 2);
pvals_ll_1y_pc = nan(size(zBeh_zCtr, 2), 2);
R2_1y_pc       = nan(size(zBeh_zCtr, 2), 3); % [reduced (static) | reduced (dyn) | unreduced]

for it = 1 : length(pvals_ll_1y_pc)
    
    Y = (zBeh_zCtr(2*Ns+1 : end, it)-zBeh_zCtr(1 : Ns, it)) ./ ...
        abs(zBeh_zCtr(1 : Ns, it));
    
    [~, p] = corr(Dyn_PCs.Scores_2w(:, 1 : nPC), Y, 'rows', 'pairwise');
    
    % Restricted models (static)
    [b, dev, stats] = glmfit(Static_info.ST, Y);
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    rL1             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_r1 = glmval(b, Static_info.ST, 'identity');
    R2_1y_pc(it, 1) = 1 - nansum((Y - yhat_r1) .^2) / nansum((Y - nanmean(Y)).^2);
    
    % Restricted models (dynamic)
    [b, dev, stats] = glmfit(Dyn_PCs.Scores_2w(:, 1 : nPC), Y);
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    rL2             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_r2 = glmval(b, Dyn_PCs.Scores_2w(:, 1 : nPC), 'identity');
    R2_1y_pc(it, 2) = 1 - nansum((Y - yhat_r2) .^2) / nansum((Y - nanmean(Y)).^2);
    
    % Unrestricted model
    [b, dev, stats] = glmfit([Dyn_PCs.Scores_2w(:, 1 : nPC), ...
        Static_info.ST], Y);
    sigma           = sqrt(nanmean(stats.resid .^ 2));
    urL             = nansum(-0.5 * (stats.resid / sigma) .^ 2 - ...
        log(sqrt(2*pi) * sigma));
    yhat_ur = glmval(b, [Dyn_PCs.Scores_2w(:, 1 : nPC), ...
        Static_info.ST], 'identity');
    R2_1y_pc(it, 3) = 1 - nansum((Y - yhat_ur) .^2) / nansum((Y - nanmean(Y)).^2);
    
    [~, pvals_ll_1y_pc(it, 1), chi2_ll_1y_pc(it, 1)] = lratiotest(urL, rL1, nPC);
    [~, pvals_ll_1y_pc(it, 2), chi2_ll_1y_pc(it, 2)] = lratiotest(urL, rL2, 1);

end

res_recovery.chi2_ll_pc = chi2_ll_1y_pc;
res_recovery.pvals_ll_pc = pvals_ll_1y_pc; 
res_recovery.R2_pc = R2_1y_pc;
res_recovery.task_list = task_list;






    

