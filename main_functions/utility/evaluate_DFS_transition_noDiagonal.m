function [DFS_trans, idx_trans, DFS_trans_v2] = evaluate_DFS_transition_noDiagonal(DFS_sw, K)

DFS_trans  = zeros(K, K);
idx_trans = [];

DFS_jumps = zeros(K, 1);

for ic = 2 : length(DFS_sw)
    k0 = DFS_sw(ic - 1);
    k1 = DFS_sw(ic);
    
    if k0 ~= k1
        DFS_jumps(k0) = DFS_jumps(k0) + 1;
        idx_trans = [idx_trans; ic];
        DFS_trans(k0, k1) = DFS_trans(k0, k1) + 1;
    end
end

DFS_trans0   = DFS_trans;
DFS_trans    = DFS_trans0 ./ sum(DFS_jumps);
DFS_trans_v2 = DFS_trans0 ./ repmat(DFS_jumps, [1 K]);
