%% GLME for fraction times
%
% inputs:
% - f_ctr: [n_CTR, K] fraction time of CTR
% - f_pat: [n_PAT, K] fraction time of PAT of one single time point (2w, 3m or 1y)
% - ST: [n_PAT, 1] 1: severe, 0: mild


function glme = GLME_FT_groups_DFS(f_ctr, f_pat, ST)

%% GLME (CTRs vs acute PAT)

disp(['GLME to test groups (CTRs and acute PATs) and DFS (1 : K) effects ', ...
    'on DFS frequency of occurrence'])

K = size(f_ctr, 2);
n_ctr = size(f_ctr, 1);
n_pat = size(f_pat, 1);

Y     = [f_ctr(:); f_pat(:)];
Yp    = round(Y * 270);

ff_ctr  = repmat(1 : K, [n_ctr, 1]);
ff_pat   = repmat(1 : K, [n_pat, 1]);
F       = categorical([ff_ctr(:); ff_pat(:)]);

g_ctr = ones(size(f_ctr));
g_2w  = repmat(ST, [1 K]);
G     = [g_ctr(:); g_2w(:)];
G(G == 1) = 2;
G(G == 0) = 3;
G = categorical(G);

s_ctr = repmat((1 : n_ctr)', [1 K]);
s_2w  = repmat(((1+n_ctr) : (n_ctr+n_pat))', [1 K]);
S     = [s_ctr(:); s_2w(:)];
S     = categorical(S);

tbl_glme = table(F, G, S, Y, Yp);

glme = fitglme(tbl_glme,'Yp ~ 1 +  F * G + (1 | S)', 'Distribution', 'Poisson');
anova(glme)