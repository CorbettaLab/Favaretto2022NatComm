%% ===================================================================== %%
%         Homotopic + DAN/DMN intra-hemispheric dynamic analysis          %
%  ===================================================================== %%

function [Homotopic_dyn_all, DAN_DMN_LR_dyn_all, Modularity_dyn_sym_neg_ave] = ...
    DFCs_homotopic_DAN_DMN_modularity(n, NET_reduction, DFC_bold, ...
    Modularity_dyn_sym_neg_ctr, Modularity_dyn_sym_neg_pat)

%% Add functions to path
addpath('utility/')

%% General Parameters

% Masks definition
mask   = triu(ones(n, n), 1);
i_mask = find(mask);

%% LOAD DATA)

Modularity_dyn_sym_neg_ctr_ave = cellfun(@(x) mean(x, 2), ...
    Modularity_dyn_sym_neg_ctr, 'UniformOutput', false);
Modularity_dyn_sym_neg_ctr_ave = cat(1, Modularity_dyn_sym_neg_ctr_ave{:});

Modularity_dyn_sym_neg_2w_ave = cellfun(@(x) mean(x, 2), ...
    Modularity_dyn_sym_neg_pat(:, 1), 'UniformOutput', false);
Modularity_dyn_sym_neg_2w_ave = cat(1, Modularity_dyn_sym_neg_2w_ave{:});

Modularity_dyn_sym_neg_3m_ave = cellfun(@(x) mean(x, 2), ...
    Modularity_dyn_sym_neg_pat(:, 2), 'UniformOutput', false);
Modularity_dyn_sym_neg_3m_ave = cat(1, Modularity_dyn_sym_neg_3m_ave{:});

Modularity_dyn_sym_neg_1y_ave = cellfun(@(x) mean(x, 2), ...
    Modularity_dyn_sym_neg_pat(:, 3), 'UniformOutput', false);
Modularity_dyn_sym_neg_1y_ave = cat(1, Modularity_dyn_sym_neg_1y_ave{:});

Modularity_dyn_sym_neg_ave = cat(1, Modularity_dyn_sym_neg_ctr_ave, ...
    Modularity_dyn_sym_neg_2w_ave, Modularity_dyn_sym_neg_3m_ave, ...
    Modularity_dyn_sym_neg_1y_ave);
Modularity_dyn_sym_neg_ave = Modularity_dyn_sym_neg_ave';

nT          = size(DFC_bold, 2);

%% HOMOTOPIC ANALYSIS (cortical w/o none)

[~, L_net, R_net, homotopic, ~] = ...
    net_mask_ind_reduced(NET_reduction, i_mask);

Homotopic_dyn_all = nan(length(homotopic.mask), nT);

for in = 1 : length(homotopic.mask)
    idx_net                     = homotopic.mask{in};
    Homotopic_dyn_all(in, :)    = nanmean(DFC_bold(idx_net, :), 1);
end

%% DAN-DMN Intrahemispheric

DAN_idx = 6;
DMN_idx = 8;

idx_DAN_DMN_L = L_net.mask{DAN_idx, DMN_idx};
idx_DAN_DMN_R = R_net.mask{DAN_idx, DMN_idx};

DAN_DMN_L_dyn_all = nanmean(DFC_bold(idx_DAN_DMN_L, :), 1);
DAN_DMN_R_dyn_all = nanmean(DFC_bold(idx_DAN_DMN_R, :), 1);

DAN_DMN_LR_dyn_all = mean([DAN_DMN_L_dyn_all; DAN_DMN_R_dyn_all], 1);

%% Remove functions to path
rmpath('utility/')
        
        
        




