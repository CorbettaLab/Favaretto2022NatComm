%% ===================================================================== %%
%           EVALUATION OF DYNAMICAL FUNCTIONAL CONNECTIVITY (DFC) AND     %
%                            LEADING EIGENVECTORS                         %
%  ===================================================================== %%
%
% inputs: 
% - n_CTR: number of control subjects
% - n_PAT: number of patients at the sub-acute stage
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
% - LEigs_bold: numerical matrix [n, nW=number of total sliding windows]:
%                 each column is the Leading Eigenvector of the DFC matrix 
%                 evaluated in the corresponding sliding window 
% - DFC_bold: numerical matrix [n*(n-1)/2, nW]: each column is the unrolled
%                 upper triangular part of the DFC matrix evaluated in the
%                 corresponding sliding window
% - Time_ctr_pat_bold: numerical matrix [3, nW]: this matrix is used to keep
%         track of the subject whose each sliding window refers to.
%     - row 1: integer indicating the group (1: CTR, 2: PAT(2w); 
%                                             3: PAT(3m); 4: PAT(1y)
%     - row 2: integer indicating the subject within each group (e.g. for a
%                 CTR subject, this value is between 1 and n_CTR
%     - row 3: integer indicating the index of the subject in the original
%                 dataset (before selection of subjects suitable for the 
%                 dynamic analysis) 
% - idx_CTR1: indexes of selected CTR from 1st session
% - idx_CTR2: indexes of selected CTR from 2nd session
% - idx_PAT_all: index of selected PAT (all time points)
% - n_Subjects: (vector [1, 4]) number ofselected subjects in each group
%                                 (CTR, PAT(2w), PAT(3m), PAT(1y))
% - idx_subj_cond: (vector [sum(n_Subjects), 1]) integer indicating the group 
%     of each subject (1: CTR, 2: PAT(2w), 3: PAT(3m), 4: PAT(1y))
%  ===================================================================== %%

function [LEs_bold, DFC_bold, Time_ctr_pat_bold, idx_CTR1, idx_CTR2, ...
    idx_PAT_all, n_Subjects, idx_subj_cond] = ...
    eval_LEs_DFC(n_CTR, n_PAT, BOLD_red_ctr1, BOLD_red_ctr2, BOLD_red_2w, BOLD_red_3m, BOLD_red_1y, ...
    NET_reduction, Tmin)

%% add functions to
addpath(genpath('utility/'))

n       = length(NET_reduction.net_reduced);        % # ROIs

n_Cond  = 3;    % 2w - 3m - 1y

% Masks definition
mask   = triu(ones(n, n), 1);
i_mask = find(mask);
l_mask = length(i_mask);

%% Select subjects with desired data

PATs_isPresent  = false(n_PAT, n_Cond);
idx_PAT         = cell(n_Cond,1 );
PATs_map        = zeros(n_PAT, n_Cond);  % map original pat to row
PATs_BOLD       = cell(n_Cond, n_PAT);

CTRs_isPresent  = false(n_CTR, 2);
CTRs_map        = zeros(n_CTR, 2);
CTRs_BOLD       = cell(2, n_CTR);

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
    
    num_pat = length(S);
    
    
    for is = 1 : num_pat
        if ~isempty(S(is).BOLD_all)
            
            BOLD_tmp  = cat(1, S(is).BOLD_all{:});
            tmask_tmp = cat(2, S(is).tmask_all{:});
            
            idx_tmask = find(tmask_tmp);
            
            if length(idx_tmask) >= Tmin
                
                pat_id = str2double(S(is).SubjectID(5 : 7));
                
                PATs_isPresent(pat_id, ic)  = true;
                PATs_map(pat_id, ic)        = is;
                
                PATs_BOLD{ic, pat_id}   = BOLD_tmp(idx_tmask(1 : Tmin), :);
            end
            
        end
    end
    
    idx_PAT{ic} = find(PATs_isPresent(:, ic));
end

idx_PAT_all         = find(sum(PATs_isPresent, 2) == 3);
BOLD_PAT_all        = PATs_BOLD(:, idx_PAT_all);
n_PAT_all           = length(idx_PAT_all);

% CONTROLS
for ic = 1 : 2
    switch ic
        case 1  
            S = BOLD_red_ctr1;
        case 2 
            S = BOLD_red_ctr2;
    end
    
    num_ctr = length(S);

    for is = 1 : num_ctr
        if ~isempty(S(is).BOLD_all)
            
            BOLD_tmp  = cat(1, S(is).BOLD_all{:});
            tmask_tmp = cat(2, S(is).tmask_all{:});
            
            idx_tmask = find(tmask_tmp);

            if length(idx_tmask) >= Tmin
                ctr_id = str2double(S(is).SubjectID(5 : 7));

                CTRs_isPresent(ctr_id, ic)  = true;
                CTRs_map(ctr_id, ic)        = is;

                CTRs_BOLD{ic, ctr_id}   = BOLD_tmp(idx_tmask(1 : Tmin), :);
            end

        end
    end
end

idx_CTR1 = find(CTRs_isPresent(:, 1));
idx_CTR2 = find(CTRs_isPresent(:, 2));
idx_CTR  = cat(1, idx_CTR1, idx_CTR2); 
n_CTR    = length(idx_CTR);

BOLD_CTR_all	= [CTRs_BOLD(1, idx_CTR1), CTRs_BOLD(2, idx_CTR2)];

%% DATA CONCATENATION

disp('=====   CONCATENATE DATA   =====')

% CONTROLS
Time_ctr  = [1 1 idx_CTR(1)]' .* ones(3, size(BOLD_CTR_all{1}, 1));

for is = 2 : n_CTR
    Time_ctr	= cat(2, Time_ctr, ...
        [1 is idx_CTR(is)]' .* ones(3, size(BOLD_CTR_all{is}, 1)));
end

BOLD_CTR    = vertcat(BOLD_CTR_all{:});

disp('CTRs ok')

% PATIENTS (2W)
Time_pat1	= [2 1 idx_PAT_all(1)]' .* ones(3, size(BOLD_PAT_all{1, 1}, 1));

for is = 2 : n_PAT_all
    Time_pat1 	= cat(2, Time_pat1, ...
        [2 is idx_PAT_all(is)]' .* ones(3, size(BOLD_PAT_all{1, is}, 1)));
end

BOLD_PATs1  = vertcat(BOLD_PAT_all{1, :});

disp('PATs 2 weeks ok')

% PATIENTS (3M)
Time_pat2	= [3 1 idx_PAT_all(1)]' .* ones(3, size(BOLD_PAT_all{2, 1}, 1));

for is = 2 : n_PAT_all
    Time_pat2 	= cat(2, Time_pat2, ...
        [3 is idx_PAT_all(is)]' .* ones(3, size(BOLD_PAT_all{2, is}, 1)));
end

BOLD_PATs2  = vertcat(BOLD_PAT_all{2, :});

disp('PATs 3 months ok')

% PATIENTS (1Y)
Time_pat3	= [4 1 idx_PAT_all(1)]' .* ones(3, size(BOLD_PAT_all{3, 1}, 1));

for is = 2 : n_PAT_all
    Time_pat3 	= cat(2, Time_pat3, ...
        [4 is idx_PAT_all(is)]' .* ones(3, size(BOLD_PAT_all{3, is}, 1)));
end

BOLD_PATs3  = vertcat(BOLD_PAT_all{3, :});

disp('PATs 1 year ok')

BOLD_PAT    = vertcat(BOLD_PATs1, BOLD_PATs2, BOLD_PATs3);
Time_pat	= cat(2, Time_pat1, Time_pat2, Time_pat3);

BOLD_ctr_pat    = vertcat(BOLD_CTR, BOLD_PAT);
Time_ctr_pat	= cat(2, Time_ctr, Time_pat);

n_Subjects = [n_CTR, n_PAT_all, n_PAT_all, n_PAT_all];

idx_subj_cond = zeros(1, sum(n_Subjects));
for ic = 1 : (n_Cond + 1)
    if ic == 1
        idx = 1 : n_Subjects(1);
    else
        idx = sum(n_Subjects(1 : (ic-1)))+1 : sum(n_Subjects(1 : ic));
    end
    
    idx_subj_cond(idx) = ic;
end

%% BOLD DFC EXTRACTION

disp('=====   BOLD DFC and LEigs EXTRACTION   =====')

win_len = 30;  % time point in each window

[Time_ctr_pat_bold, LEs_bold, ~, DFC_bold] = ...
    bold_DFC_LEs_extraction(n_Subjects, n_Cond, Time_ctr_pat, n, Tmin, win_len, ...
    i_mask, l_mask, BOLD_ctr_pat);

%% remove functions from path
rmpath(genpath('utility/'))



