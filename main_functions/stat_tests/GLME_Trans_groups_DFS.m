%% GLME for transitions (CTR vs PAT)
% 
% inputs:
% - T_ctr_vec: [n*(n-1), n_CTR] vectorized Trans matrix of CTR
% - T_pat_vec: [n*(n-1), n_PAT] vectorized Trans matrix of PAT (one single time point)
% - ST [n_PAT, 1] 1: severe, 0: mild
% - K: number of DFSs


function glme_t = GLME_Trans_groups_DFS(T_ctr_vec, T_pat_vec, ST, K)

mask_K      = ones(K, K) - eye(K, K);
mask_K_ind  = find(mask_K);

n_ctr = size(T_ctr_vec, 2);
n_pat = size(T_pat_vec, 2);

Y = [T_ctr_vec, T_pat_vec]';
Y = Y(:);

TT = repmat((1 : length(mask_K_ind)), [n_ctr+n_pat, 1]);
TT = TT(:);
TT = categorical(TT);

G_pat = nan(n_pat, 1);
G_pat(ST == 1) = 2;
G_pat(ST == 0) = 3;
G = [ones(n_ctr, 1); G_pat];

G = repmat(G, [1 length(mask_K_ind)]);
G = G(:);
G = categorical(G);

S = repmat((1 : (n_ctr+n_pat))', [1 length(mask_K_ind)]);
S = S(:);
S = categorical(S);

tbl_gmle_t = table(TT, G, S, Y);
glme_t = fitglme(tbl_gmle_t, 'Y ~ 1 + TT * G + (1 | S)');

anova(glme_t)