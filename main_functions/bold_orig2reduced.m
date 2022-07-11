%% ===================================================================== %%
%       TIME SERIES DATA TRANSFORMATION FROM ORIGINAL TO REDUCED SPACE    %
%  ===================================================================== %%
% 
% Inputs:
% - NET_reduction (output of "dimensionality_reduction.m")
% - BOLD_runs: (Structure [n_subj, 1]) original BOLD data.
%     For each subject:
%         - SubjectID: string with subject ID
%         - BOLD_all: (Cell [n_run, 1]) BOLD data for each run
%         - tmask_all: (Cell [1, n_run]) boolean vector indicating good/bad
%         frames in each run
% 
% Output:
% - BOLD_reduced_runs: as input 'BOLD_runs' in the reduced space with an
% additional field for each subject:
%     - sFC_all: static FC evaluated for each run


function BOLD_reduced_runs = bold_orig2reduced(NET_reduction, BOLD_runs)

%% add functions to path
addpath('utility/')

%% Reduce dimensionality

n_SUB = length(BOLD_runs);

BOLD_reduced_runs = struct('SubjectID', cell(1, n_SUB), ...
    'BOLD_all', cell(1, n_SUB), 'sFC_all', cell(1, n_SUB), ...
    'tmask_all', cell(1, n_SUB));

for is = 1 : n_SUB
    
    BOLD_reduced_runs(is).SubjectID = BOLD_runs(is).SubjectID;
    BOLD_reduced_runs(is).tmask_all = BOLD_runs(is).tmask_all;
    
    
    BOLD_reduced_runs(is).BOLD_all = ...
        reduced_Runs(BOLD_runs(is).BOLD_all, NET_reduction.final_index);
    
    n_runs = length(BOLD_reduced_runs(is).BOLD_all);
    
    BOLD_reduced_runs(is).sFC_all = cell(n_runs, 1);
    
    
    
    for ir = 1 : n_runs
        idx_tmask = find(BOLD_reduced_runs(is).tmask_all{ir});
        
        if length(idx_tmask) > 10
            BOLD_reduced_runs(is).sFC_all{ir} = ...
                corr(BOLD_reduced_runs(is).BOLD_all{ir}(idx_tmask, :));
        end
    end
    
    disp(['SUB ', num2str(is), ': ok!'])
    
    
end

%% remove functions to path
rmpath('utility/')












