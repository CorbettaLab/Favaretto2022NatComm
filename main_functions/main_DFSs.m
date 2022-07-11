%% ===================================================================== %%
%                        DYNAMICAL FUNCTIONAL STATES                      %
%  ===================================================================== %%
% 
% inputs:
% - LEs_bold: (matrix: [n, nW]) (output of "eval_LEs_DFC.m")
% - K: (integer >= 2) number of DFSs
% 
% outputs:
% - DFS_sw: (vector: [nW, 1]) vector indicating to which DFS each
%                              sliding window belongs
% - DFS: (matrix: [n*(n-1)/2, K]) unrolled upper triangular part each DFS matrix
% - DFS_dist: (matrix: [nW, K]) distance from each DFS of each sliding window
% - DFC_LEs: (matrix [n*(n-1)/2, nW]): each column is the unrolled
%                 upper triangular part of the matrix evaluated as V*V',
%                 where V is the Leading eigenvector evaluated at the related
%                 sliding window
%  ===================================================================== %%


function [DFS_sw, DFS, DFS_dist, DFC_LEs] = main_DFSs(LEs_bold, K)

%% Add toolbox to path

addpath(genpath('utility/'))

%% LEiDA (linkage, ward)

% create DFC matrix based only on the first eigenvector: DFC_LE = V * V'
[n, nT] = size(LEs_bold);

% Masks definition
mask   = triu(ones(n, n), 1);
i_mask = find(mask);
l_mask = length(i_mask);

DFC_LEs = nan(l_mask, nT);

for it = 1 : nT
    DFC_it              = LEs_bold(:, it) * (LEs_bold(:, it)');
    DFC_LEs(:, it)    = DFC_it(i_mask);
end

[DFS_sw, DFS, ~, DFS_dist] = kmeans(DFC_LEs', K, 'Replicates', 20, ...
        'MaxIter', 200, 'Distance', 'correlation');

DFS = DFS';
  

%% Remove toolbox from path

rmpath(genpath('utility/'))

