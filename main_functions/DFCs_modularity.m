%% ===================================================================== %%
%                       Modularity dynamic analysis                       %
%  ===================================================================== %%

function [Modularity_dyn_sym_neg_ctr, Modularity_dyn_sym_neg_pat] = ...
    DFCs_modularity(n, NET_reduction, DFC_bold, n_CTR, n_PAT, Time_ctr_pat_bold)

%% Add functions to path
addpath(genpath('utility'))
addpath('external/BCT')

%% General Parameters

% Masks definition
mask   = triu(ones(n, n), 1);
i_mask = find(mask);
l_mask = length(i_mask);

[r_i, c_i] = ind2sub([n n], i_mask);


%% Modularity (with a priori modules)
%  Consider only cortical ROIs (no NONE networks)

idx_mod = find(NET_reduction.net_reduced <= 8);

net = NET_reduction.net_reduced(idx_mod);

% mask
mask_mod   = triu(ones(length(idx_mod)), 1);
i_mask_mod = find(mask_mod);
l_mask_mod = length(i_mask_mod);

% Density definition
prc     = 4 : 4 : 50;
n_prc   = length(prc);
val_prc = ceil(l_mask_mod .* .5 .* prc ./ 100);


%% Symmetric Negative modularity

Modularity_dyn_sym_neg_ctr = cell(n_CTR, 1);
Modularity_dyn_sym_neg_pat = cell(n_PAT, 3);

for is = 1 : n_CTR
    
    idx_is = find((Time_ctr_pat_bold(1, :) == 1) & ...
        (Time_ctr_pat_bold(2, :) == is));
    DFC_is = DFC_bold(:, idx_is)';
    
    Modularity_dyn_sym_neg_ctr{is} = nan(length(idx_is), n_prc);
    
    for it = 1 : length(idx_is)
        DFC_vec = DFC_is(it, :)';
        mat = zeros(n, n);
        for ii = 1 : l_mask
            mat(r_i(ii), c_i(ii)) = DFC_vec(ii);
        end
        mat = mat + mat';
        mat = mat(idx_mod, idx_mod);
        mat(isinf(mat)) = 0;
        mat(isnan(mat)) = 0;
        
        val_sort = sort(abs(mat(i_mask_mod)), 'descend');
        
        for ip = 1 : n_prc
            val_ip = val_sort(val_prc(ip));
            mat_ip = mat;
            mat_ip(abs(mat_ip) < val_ip) = 0;
            
            Modularity_dyn_sym_neg_ctr{is}(it, ip) = ...
                community_louvain_apriori(mat_ip, 1, net, 'negative_sym');
        end
    end
end

for cond = 1 : 3
    for is = 1 : n_PAT
        
        idx_is = find((Time_ctr_pat_bold(1, :) == (cond+1)) & ...
            (Time_ctr_pat_bold(2, :) == is));
        DFC_is = DFC_bold(:, idx_is)';
        
        Modularity_dyn_sym_neg_pat{is, cond} = nan(length(idx_is), n_prc);
        
        for it = 1 : length(idx_is)
            DFC_vec = DFC_is(it, :)';
            mat = zeros(n, n);
            for ii = 1 : l_mask
                mat(r_i(ii), c_i(ii)) = DFC_vec(ii);
            end
            mat = mat + mat';
            mat = mat(idx_mod, idx_mod);
            mat(isinf(mat)) = 0;
            mat(isnan(mat)) = 0;
            
            val_sort = sort(abs(mat(i_mask_mod)), 'descend');
            
            for ip = 1 : n_prc
                val_ip = val_sort(val_prc(ip));
                mat_ip = mat;
                mat_ip(abs(mat_ip) < val_ip) = 0;
                
                Modularity_dyn_sym_neg_pat{is, cond}(it, ip) = ...
                    community_louvain_apriori(mat_ip, 1, net, 'negative_sym');
            end
        end
    end
end

%% Remove functions from path
rmpath(genpath('utility'))
rmpath('external/BCT')
