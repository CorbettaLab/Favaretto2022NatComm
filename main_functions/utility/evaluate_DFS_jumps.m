function [DFS_trans, idx_trans] = evaluate_DFS_jumps(DFS_sw, K)

DFS_trans  = zeros(K, K);
idx_trans = [];

for ic = 2 : length(DFS_sw)
    k0 = DFS_sw(ic - 1);
    k1 = DFS_sw(ic);
    
    if k0 ~= k1
        idx_trans = [idx_trans; ic];
        DFS_trans(k0, k1) = DFS_trans(k0, k1) + 1;
    end
end
