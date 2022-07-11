%% ===================================================================== %%
%                          create_data_example.m                          %
%  Important Notes:                                                       %
%  use this script to generate an example dataset that contains random    %
%  generated time series for two subjects. The aim of this dataset is     %
%  ONLY to show the formal structure required for input data.             %
%  This dataset DO NOT simulate real BOLD fMRI time series.               %
%  (Data are randomly generated without any assumptions)                  %
%                                                                         %
%  output:                                                                %
%  - DATA: structure of size [1, n_subjects] and 3 fields:                %
%           1. SubjectID: string indicating the anonymous ID of each      %
%           subject (es. 'FCS_001_AMC')                                   %
%           2. BOLD_all: cell of size [n_runs, 1].                        %
%                        Each run is matrix of size [T, n_ROIs]           %
%           3. tmask_all: cell of size [n_runs, 1].                       %
%                         Each run is boolean vector of size [1, T],      %
%                         indicating if a frame is of good (1) or bad (0) %
%                         quality                                         %
%                                                                         %
%  ===================================================================== %%

clc
close all
clear all

T_run   = 128; % TR in each run
n       = 343; % number of ROIs
max_bad = 20; % max number of bad frame for each run

DATA = struct('SubjectID', {'FCS_001_AMC', 'FCS_002_AMC'}, ...
    'BOLD_all', {cell(7, 1), cell(7, 1)}, ...
    'tmask_all', {cell(7, 1), cell(7, 1)});

for ii = 1 : 2
    r = length(DATA(ii).BOLD_all);
    for ir = 1 : r
        % random BOLD data
        DATA(ii).BOLD_all{ir} = rand(T_run, n);
        % random tmask
        DATA(ii).tmask_all{ir} = ones(1, T_run);
        
        idx_perm = randperm(T_run);
        n_bad = randi(max_bad);
        DATA(ii).tmask_all{ir}(idx_perm(1 : n_bad)) = 0;
    end
end

save DATA_example.mat DATA
        