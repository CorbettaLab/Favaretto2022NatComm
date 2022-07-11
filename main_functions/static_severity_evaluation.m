%% ===================================================================== %%
%     Static Principal Component (ST) evaluation to split patients in     %
%                SEVERE and MILD based on static FC impaiment             %
%  ===================================================================== %%
% 
% inputs: 
% - n_CTR: number of control subjects
% - n_PAT: number of patients at the sub-acute stage
% - idx_PAT_all: index of selected PAT (all time points) (output of
% "eval_LEs_DFC.m")
% - BOLD_red_ctr1: (Struct) (output of "bold_orig2reduced.m" for 
%                                             1st session of CTRs subjects)
% - BOLD_red_ctr2: (Struct) (output of "bold_orig2reduced.m" for 
%                                             2nd session of CTRs subjects)
% - BOLD_red_2w: (Struct) (output of "bold_orig2reduced.m" for PATs at 2 weeks)
% - BOLD_red_3m: (Struct) (output of "bold_orig2reduced.m" for PATs at 3 months)
% - BOLD_red_1y: (Struct) (output of "bold_orig2reduced.m" for PATs at 1 year)
% - NET_reduction (Struct) (output od "dimensionality_reduction.m")
% - Tmin: minimum number of good frames required to evaluate the Dynamical FC
% 
% outputs:
% - Static_info: (Struct) all information related to static analysis, in
% particular Static_info.ST contains is a vector [n_PAT_all, 1] such that:
%     - Static_info.ST(i) >= 0 => severe
%     - Static_info.ST(i) < 0 => mild
%  ===================================================================== %%
                   

function Static_info = static_severity_evaluation(n_CTR, n_PAT, idx_PAT_all, ...
    BOLD_red_ctr1, BOLD_red_ctr2, BOLD_red_2w, BOLD_red_3m, BOLD_red_1y, ...
    NET_reduction, Tmin)


%% Add toolbox to path

addpath(genpath('utility/'))
addpath(genpath('external/NIfTI_20140122/'))
addpath('external/BCT')

%% General Parameters

n       = length(NET_reduction.net_reduced);        % # ROIs

n_Cond  = 3;    % 2w - 3m - 1y

% Masks definition
mask   = triu(ones(n, n), 1);
i_mask = find(mask);
l_mask = length(i_mask);


PATs_isPresent  = false(n_PAT, n_Cond);
idx_PAT         = cell(n_Cond,1 );
PATs_map        = zeros(n_PAT, n_Cond);  % map original pat to row
PATs_sFC        = cell(n_Cond, 1);

CTRs_isPresent  = false(n_CTR, 2);
CTRs_map        = zeros(n_CTR, 2);
CTRs_sFC        = cell(2, 1);

% PATIENTS
for ic = 1 : n_Cond
    switch ic
        case 1  % 2 weeks
            S = BOLD_red_2w;
        case 2  % 3 months
            S = BOLD_red_3m;
        case 3  % 1 year
            S = BOLD_red_1y;
    end
    
    num_pat     = length(S);
    
    PATs_sFC{ic} = nan(n, n, n_PAT);
    
    for is = 1 : num_pat

        sFC_is  = S(is).sFC_all;
        l       = cellfun(@(x) length(x), sFC_is);
        tmask = cat(2, S(is).tmask_all{:});
        good_frames = length(find(tmask == 1));
        
        if (~isempty(find(l > 0))) && (good_frames >= Tmin)
            
            sFC_tmp = S(is).sFC_all(find(l > 0));
            sFC_tmp = cat(3, sFC_tmp{:});
            
            pat_id = str2double(S(is).SubjectID(5 : 7));
            
            PATs_sFC{ic}(:, :, pat_id) = nanmean(sFC_tmp, 3);

            PATs_isPresent(pat_id, ic)  = true;
            PATs_map(pat_id, ic)        = is;
        end
    end
    
    idx_PAT{ic} = find(PATs_isPresent(:, ic));
end

n_PAT = sum(PATs_isPresent, 1); % # of patients at each time point

PATs_sFC{1}        = PATs_sFC{1}(:, :, idx_PAT{1});
PATs_sFC{2}        = PATs_sFC{2}(:, :, idx_PAT{2});
PATs_sFC{3}        = PATs_sFC{3}(:, :, idx_PAT{3});


% CONTROLS
for ic = 1 : 2
    switch ic
        case 1  
            S = BOLD_red_ctr1;
        case 2 
            S = BOLD_red_ctr2;
    end
    
    num_ctr                 = length(S);

    CTRs_sFC{ic} = nan(n, n, n_CTR);
    
    for is = 1 : num_ctr
        sFC_is  = S(is).sFC_all;
        l       = cellfun(@(x) length(x), sFC_is);
        
        tmask = cat(2, S(is).tmask_all{:});
        good_frames = length(find(tmask == 1));
        
        if (~isempty(find(l > 0))) && (good_frames >= Tmin)
            
            sFC_tmp = S(is).sFC_all(find(l > 0));
            sFC_tmp = cat(3, sFC_tmp{:});
            
            ctr_id = str2double(S(is).SubjectID(5 : 7));

            CTRs_sFC{ic}(:, :, ctr_id) = nanmean(sFC_tmp, 3);
            
            CTRs_isPresent(ctr_id, ic)  = true;
            CTRs_map(ctr_id, ic)        = is;

        end
    end
end

idx_CTR1 = find(CTRs_isPresent(:, 1));
idx_CTR2 = find(CTRs_isPresent(:, 2));
idx_CTR  = cat(1, idx_CTR1, idx_CTR2); 
n_CTR   = length(idx_CTR);

CTRs_sFC = cat(3, CTRs_sFC{1}(:, :, idx_CTR1), CTRs_sFC{2}(:, :, idx_CTR2));


% z-Fisher Transform
zFC_ctr    = atanh(CTRs_sFC);
zFC_pat_2w = atanh(PATs_sFC{1});
zFC_pat_3m = atanh(PATs_sFC{2});
zFC_pat_1y = atanh(PATs_sFC{3});

%% FC vectors

FC_vec_ctr   = nan(l_mask, n_CTR);

for is = 1 : n_CTR
    FC_tmp = squeeze(zFC_ctr(:, :, is));
    FC_tmp(isinf(FC_tmp)) = 0;
    
    FC_vec_ctr(:, is) = FC_tmp(i_mask);   
end

FC_vec_2w   = nan(l_mask, n_PAT(1));

for is = 1 : n_PAT(1)
    FC_tmp = squeeze(zFC_pat_2w(:, :, is));
    FC_tmp(isinf(FC_tmp)) = 0;
    FC_tmp(isnan(FC_tmp)) = 0;
    
    FC_vec_2w(:, is) = FC_tmp(i_mask);
end

FC_vec_3m   = nan(l_mask, n_PAT(2));

for is = 1 : n_PAT(2)
    FC_tmp = squeeze(zFC_pat_3m(:, :, is));
    FC_tmp(isinf(FC_tmp)) = 0;
    FC_tmp(isnan(FC_tmp)) = 0;
    
    FC_vec_3m(:, is) = FC_tmp(i_mask);
end

FC_vec_1y   = nan(l_mask, n_PAT(3));

for is = 1 : n_PAT(3)
    FC_tmp = squeeze(zFC_pat_1y(:, :, is));
    FC_tmp(isinf(FC_tmp)) = 0;
    FC_tmp(isnan(FC_tmp)) = 0;
    
    FC_vec_1y(:, is) = FC_tmp(i_mask); 
end

FC_vec_all      = cat(2, FC_vec_ctr, FC_vec_2w, FC_vec_3m, FC_vec_1y);

%% Z-Scored w.r.t. CTR

zFC_ctr_ave = nanmean(zFC_ctr, 3);
zFC_ctr_std = nanstd(zFC_ctr, [], 3);

FC_pat_2w_Zctr = (zFC_pat_2w - repmat(zFC_ctr_ave, [1 1 n_PAT(1)])) ./ ...
    repmat(zFC_ctr_std, [1 1 n_PAT(1)]);

FC_pat_2w_Zctr_vec = nan(l_mask, n_PAT(1));

for is = 1 : n_PAT(1)
    data = squeeze(FC_pat_2w_Zctr(:, :, is));
    FC_pat_2w_Zctr_vec(:, is) = data(i_mask);
end

%% PCA of Z-Scored values

[~, idx_noDyn] = setdiff(idx_PAT{1}, idx_PAT_all);
[~, idx_Dyn]   = intersect(idx_PAT{1}, idx_PAT_all);


%% PCA of Z-Scored values

[PC_2w_ctr, score_2w_ctr_noDyn, Eig2_2w_ctr, ~, expVar_2w_ctr, mu_2w_ctr] = ...
    pca(FC_pat_2w_Zctr_vec(:, idx_noDyn)');

% Change sign of PC1 and PC2
PC_2w_ctr(:, 2) = - PC_2w_ctr(:, 2);
score_2w_ctr_noDyn(:, 2) = - score_2w_ctr_noDyn(:, 2);

% Reproject Dynamical PATs on the same space
data = FC_pat_2w_Zctr_vec(:, idx_Dyn)';
score_2w_ctr = nan(length(idx_PAT_all), size(score_2w_ctr_noDyn, 2));

for ii = 1 : length(idx_PAT_all)
    X = data(ii, :) - mu_2w_ctr;
    idx_nn = find(~isnan(X));
    score_2w_ctr(ii, :) = X(idx_nn) * PC_2w_ctr(idx_nn, :);
end

%% ST definition
ST = nansum(score_2w_ctr(:, 1 : 2), 2);

% PCs matrices
[ri, ci] = ind2sub([n n], i_mask);

PCs_2w_ctr_mat = zeros(n, n, size(PC_2w_ctr, 2));

for ii = 1 : l_mask
    PCs_2w_ctr_mat(ri(ii), ci(ii), :) = PC_2w_ctr(ii, :);
    PCs_2w_ctr_mat(ci(ii), ri(ii), :) = PC_2w_ctr(ii, :);
end


%% Homotopic (cortical w/o none) connections

idx_sbj_cond = cat(1, ones(n_CTR, 1), 2*ones(n_PAT(1), 1), ... 
    3*ones(n_PAT(2), 1), 4*ones(n_PAT(3), 1));


[~, L_net, R_net, homotopic, ~] = ...
    net_mask_ind_reduced(NET_reduction, i_mask);

Homotopic_all = nan(length(homotopic.mask), length(idx_sbj_cond));

for in = 1 : length(homotopic.mask)
    idx_net                 = homotopic.mask{in};
    Homotopic_all(in, :)    = nanmean(FC_vec_all(idx_net, :), 1);
end

%% DAN-DMN Intrahemispheric

DAN_idx = 6;
DMN_idx = 8;

idx_DAN_DMN_L = L_net.mask{DAN_idx, DMN_idx};
idx_DAN_DMN_R = R_net.mask{DAN_idx, DMN_idx};

DAN_DMN_L_all = nanmean(FC_vec_all(idx_DAN_DMN_L, :), 1);
DAN_DMN_R_all = nanmean(FC_vec_all(idx_DAN_DMN_R, :), 1);

    
%% Modularity (only cortical - no none)

idx_mod = find(NET_reduction.net_reduced <= 8);
net     = NET_reduction.net_reduced(idx_mod);

% mask
mask_mod   = triu(ones(length(idx_mod)), 1);
i_mask_mod = find(mask_mod);
l_mask_mod = length(i_mask_mod);

% Density definition
prc     = 4 : 4 : 50;
n_prc   = length(prc);
val_prc = ceil(l_mask_mod .* .5 .* prc ./ 100);


% Only positive modularity
Modularity_static_pos_ctr = nan(n_CTR, n_prc);
Modularity_static_pos_pat = cell(3, 1);

for is = 1 : n_CTR
    mat = squeeze(zFC_ctr(:, :, is));
    mat = mat(idx_mod, idx_mod);
    mat(isinf(mat)) = 0;
    mat(isnan(mat)) = 0;
    mat(mat < 0)    = 0;
    
    val_sort = sort(mat(i_mask_mod), 'descend');
    
    for ip = 1 : n_prc
        val_ip = val_sort(val_prc(ip));
        mat_ip = mat;
        mat_ip(mat_ip < val_ip) = 0;
        
        Modularity_static_pos_ctr(is, ip) = ...
            community_louvain_apriori(mat_ip, 1, net, 'modularity');
    end
end

for cond = 1 : 3
    
    Modularity_static_pos_pat{cond} = nan(n_PAT(cond), n_prc);
    
    switch cond
        case 1
            FC_pat = zFC_pat_2w;
        case 2
            FC_pat = zFC_pat_3m;
        case 3
            FC_pat = zFC_pat_1y;
    end
    for is = 1 : n_PAT(cond)
        mat = squeeze(FC_pat(:, :, is));
        mat = mat(idx_mod, idx_mod);
        mat(isinf(mat)) = 0;
        mat(isnan(mat)) = 0;
        mat(mat < 0)    = 0;
        
        val_sort = sort(mat(i_mask_mod), 'descend');
        
        for ip = 1 : n_prc
            val_ip = val_sort(val_prc(ip));
            mat_ip = mat;
            mat_ip(mat_ip < val_ip) = 0;
            
            Modularity_static_pos_pat{cond}(is, ip) = ...
                community_louvain_apriori(mat_ip, 1, net, 'modularity');
            
        end
    end
end


%% Save data

Static_info.ST = ST;
Static_info.DAN_DMN_L_all = DAN_DMN_L_all;
Static_info.DAN_DMN_R_all = DAN_DMN_R_all;
Static_info.Eig2_2w_ctr = Eig2_2w_ctr; 
Static_info.FC_vec_all = FC_vec_all;
Static_info.Homotopic_all = Homotopic_all;
Static_info.PC_2w_ctr = PC_2w_ctr; 
Static_info.PCs_2w_ctr_mat = PCs_2w_ctr_mat; 
Static_info.expVar_2w_ctr = expVar_2w_ctr;
Static_info.mu_2w_ctr = mu_2w_ctr; 
Static_info.score_2w_ctr_noDyn = score_2w_ctr_noDyn; 
Static_info.score_2w_ctr = score_2w_ctr;
Static_info.idx_CTR1 = idx_CTR1; 
Static_info.idx_CTR2 = idx_CTR2;
Static_info.idx_Dyn = idx_Dyn;
Static_info.idx_noDyn = idx_noDyn;
Static_info.idx_PAT = idx_PAT;
Static_info.idx_PAT_all = idx_PAT_all;
Static_info.idx_sbj_cond = idx_sbj_cond;
Static_info.zFC_ctr = zFC_ctr;
Static_info.zFC_pat_2w = zFC_pat_2w;
Static_info.zFC_pat_3m = zFC_pat_3m;
Static_info.zFC_pat_1y = zFC_pat_1y;

%% Remove functions and toolboxes from path

rmpath(genpath('utility/'))

rmpath(genpath('external/NIfTI_20140122/'))
rmpath('external/BCT')

