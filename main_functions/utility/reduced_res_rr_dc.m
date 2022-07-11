function [beta_matrix, b_red, R2, pval] = reduced_res_rr_dc(output, ...
    mask, NET_info)



%% Define mask indexes of disconnections

n = size(mask, 1);  % # nodes
imask = find(mask);

[ri, ci] = ind2sub([n n], imask);

%% Plot Results RR

X           = output.X;
R2          = output.RR2;


% Model significance
R2_perm = output.RR2_perm;

pval = length(find(R2_perm > R2)) / length(R2_perm);



%% Plot Beta into the brain

beta_opt  = output.beta_ave;
beta_perm = output.beta_perm;

prc         = prctile(beta_perm, [2.5 97.5], 2);
idx_beta    = find(beta_opt < prc(:, 1) | beta_opt > prc(:, 2));
% idx_beta = 1 : length(beta_opt);

coeff   = output.coeff(:, 1 : output.nPC);
mPC     = output.mX;
sPC     = output.sX;

% beta_voxel_pc = beta_opt(idx_beta)' * coeff(:, idx_beta)';

beta_dcm_pc = ...
    beta_reprojection(beta_opt(idx_beta), X, coeff(:, idx_beta), mPC, sPC);

% Renormalize
cbeta_dcm_pc = beta_dcm_pc ./ max(abs(beta_dcm_pc(:)));


% Come back to the full set of disconnections
beta_dcm = nan(size(output.Xorig, 2), 1);
beta_dcm(output.idx_voxels) = beta_dcm_pc;


%% extract betas into lesion space and save

% Threshold on beta_voxel;
beta_thr = 0;
beta_dcm(abs(beta_dcm) < beta_thr) = nan;

% Rescale again
beta_dcm = beta_dcm ./ max(abs(beta_dcm(:)));

beta_matrix = nan(n, n);

for ii = 1 : length(imask)
    beta_matrix(ri(ii), ci(ii)) = beta_dcm(ii);
    beta_matrix(ci(ii), ri(ii)) = beta_dcm(ii);
end

%% Project onto the reduced parcellation

n_red = 90;

subcort_index_359 = {[9 11 13 15 17 19 21 23 25], 7, 1, 3, 5, 35, [], [], [], [], ...
    [10 12 14 16 18 20 22 24 26], 8, 2, 4, 6, [], [], [], []}';

subcort_index_359 = cellfun(@(x)  x+324, subcort_index_359, ...
    'UniformOutput', false);

final_index = cat(1, NET_info.final_index(1 : 71), subcort_index_359);

b_red = nan(n_red, n_red);

for ii = 1 : n_red
    for iii = ii : n_red
        b_tmp = beta_matrix(final_index{ii}, final_index{iii});
        b_red(ii, iii) = nanmean(b_tmp(:));
        if ii ~= iii
            b_red(iii, ii) = b_red(ii, iii);
        end
    end
end

end






