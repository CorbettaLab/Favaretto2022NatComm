function Runs_red_TS = reduced_Runs(Runs_TS, index)

Runs_red_TS = cell(size(Runs_TS));

n = length(index);

n_runs = length(Runs_TS);

for ir = 1 : n_runs
    Runs_red_TS{ir} = nan(size(Runs_TS{ir}, 1), n);
    
    for in = 1 : n
        Runs_red_TS{ir}(:, in) = nanmean(Runs_TS{ir}(:, index{in}), 2);
    end
end