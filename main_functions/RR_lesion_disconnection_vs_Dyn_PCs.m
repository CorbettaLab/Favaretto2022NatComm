%% ===================================================================== %%
%        Relationship between lesions/disconnections and Dynamical PCs    %                         
%  ===================================================================== %%
% 
% inputs:
% - Dyn_PCs: all information of Dynamical PCs (output of
% "dynamic_measures_PCA.m")
% - folder_les: folder containing the following data:
%     - 'subject_list.mat': list of subjects with lesion informatio 
% 	- 'lesion_matrix_MNI_333_fnov.mat': binaty lesion volume masks
% 	- 'parcel_discon_359_fnov.mat': number of disconnected streamlines between each pair of ROIs      
% 	- 'sc_template.mat': number of streamlines that connect each pair of ROIs in the template
% - idx_PAT_all: index of selected PAT (all time points)
% - NET_reduction: output of "dimensionality_reduction.m"
% 
% outputs:
% - rr_les: (Struct) results of Ridge Regression using lesion information to
% explain Dyn-PCs
% - rr_dcm: (Struct) results of Ridge Regression using disconnectome information to
% explain Dyn-PCs
%  ===================================================================== %%


function [rr_les, rr_dcm] = RR_lesion_disconnection_vs_Dyn_PCs(Dyn_PCs, folder_les, idx_PAT_all, NET_reduction)

%% add functions to path
addpath('utility/')


% lesion data
load(fullfile(folder_les, 'subject_list.mat'), 'subject_list')  % list of subjects with lesion informatio 
load(fullfile(folder_les, 'lesion_matrix_MNI_333_fnov.mat'), 'lesion_matrix_MNI_333_fnov')    % binaty lesion volume masks
load(fullfile(folder_les, 'parcel_discon_359_fnov.mat'), 'parcel_discon_359_fnov')        % number of disconnected streamlines between each pair of ROIs      
load(fullfile(folder_les, 'sc_template.mat'), 'sc_template')                   % number of streamlines that connect each pair of ROIs in the template

%% Select Subjects used in the dynamical states

[~, idx_les, idx_dyn] = intersect(subject_list, idx_PAT_all);

lesions          = lesion_matrix_MNI_333_fnov(idx_les, :, :, :);
disconnections   = parcel_discon_359_fnov(idx_les, :, :);

% vectorized lesions and disconnections
[Ns, d1, d2, d3] = size(lesions);
l_les = d1 * d2 * d3;

d4      = size(disconnections, 2);
mask    = triu(ones(d4, d4), 1);
i_mask  = find(mask);

SC_vec      = sc_template(i_mask);   % vectorized CS
ind_zero    = find(SC_vec > 0);

X_les = nan(Ns, l_les);
X_dcs = nan(Ns, length(i_mask));

for is = 1 : Ns
    % lesions
    les_is          = squeeze(lesions(is, :, :, :));
    X_les(is, :)    = les_is(:);
    
    % disconenctions
    dcs_is              = squeeze(disconnections(is, :, :));
    dcs_is              = dcs_is(i_mask);
    dcs_is(ind_zero)    = dcs_is(ind_zero) ./ SC_vec(ind_zero);
    X_dcs(is, :)        = dcs_is;
end

%% Run ridge regression analysis -- Lesions vs DynPCs

options_rr = getOptionsdefault;
options_rr.thr_les_size(2) = size(X_les, 2);

rr_les = cell(3, 1);

for ip = 1 : 3 % PC1, PC2, PC3

    rr_les{ip} = ridge_regression_final(X_les, Dyn_PCs.Scores_2w(:, ip), ...
        ['rr_pc', num2str(ip), '_les'], options_rr);
    
    % Plot RR results
    rr_les{ip}.name      = ['rr_pc', num2str(ip), '_les'];
    rr_les{ip}.pc        = ip;
    rr_les{ip}.R2_pval   = length(find(rr_les{ip}.RR2_perm > rr_les{ip}.RR2)) / length(rr_les{ip}.RR2_perm);

end

%% Run ridge regression analysis -- Disconection vs DynPCs

options_rr = getOptionsdefault;
options_rr.thr_les_size(2) = max(sum(X_dcs, 2));


rr_dcm = cell(3, 1);

for ip = 1 : 3 % PC1, PC2, PC3

    rr_dcm{ip} = ridge_regression_final(X_dcs, Dyn_PCs.Scores_2w(:, ip), ...
        ['rr_pc', num2str(ip), '_dcs'], options_rr);
        
    % Plot RR results
    rr_dcm{ip}.name      = ['rr_pc', num2str(ip), '_dcm'];
    rr_dcm{ip}.pc        = ip;
    
    [rr_dcm{ip}.beta_matrix, rr_dcm{ip}.beta_reduced, rr_dcm{ip}.R2_pval] = ...
        reduced_res_rr_dc(rr_dcm{ip}, mask, NET_reduction);

end


%% remove functions to path
rmpath('utility/')