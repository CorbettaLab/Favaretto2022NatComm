function lifespan = evaluate_dwelltime(DFS_sw, K)

lifespan = repmat({[]}, [K, 1]);

iter   = 0;
ic_old = 0;



for ii = 1 : length(DFS_sw)
    ic = DFS_sw(ii);
    
    if ic == ic_old
        iter = iter + 1;
    elseif ic_old > 0
        lifespan{ic_old}(end+1) = iter;
        iter = 1;
    else
        iter = 1;
    end
    
    ic_old = ic;
end
  
lifespan{ic_old}(end+1) = iter;