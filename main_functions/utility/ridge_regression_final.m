%% ===================================================================== %%
%                             ridge_regression.m                             %
%
% This function runs a principal component ridge regression analysis using
% structural lesion data to predict an outcome measure.
%
% X data: z-normalized w.r.t. the whola matrix
% Y data: centered or z-normalized
% Ridge function: without intercept
% X data preprocessing: only PCA
% number of PCs: two thresholds (95% and 99%)
% Regularization: lambda (-1e-5 : 1e+5) selected through leave-one-out
%                 always check if it exists a minimum
% Model Fit: use averaged beta (across subjects) obtained with leave-one-out
%            at optimal lambda
% Post-processing beta: z-normalization of betas before comung back to the
%                       voxels space
% Map visualization: use a threshold of 10% of the maximum value too
%                    visualize betas in the voxels space
% Model test: permutation test (10,000 permutation of behavioral scores)

% Input - X need to be a 2D matrix [Nsubjects X #voxels]
%         Y need to be a column vector [Nsubjects X 1]
%         test_name - is the name of the test that will be used to save the betas correspondent
%         output_folder - is the folder where the betas are stored in nifti_form and
%                         the workspace is saved to perform further analysis (if requested)
%         options   - is a structure with parameter to tune the ridge regression
%                   the field of the structure are reported here with the default value
%                   assigned if no option structure is provided
%                   options.nPC_thresh = 95
%                   options.l = logspace(-5, +5, 200);
%                   options.Nperm = 10000;
%                   options.save_workspace = false;
%                   options.template = template of the nifti_file (MNI_2006_1mm is the default)
%                   options.zscoreY = true
%                   options.thr_sbj_overlap = 1;


function res_all = ridge_regression_final(X, Y, task_name, options)

disp('==================================================================')
disp(['                             ', task_name])
disp('==================================================================')

if nargin < 4
    options=getOptionsdefault();
end

if isempty(task_name)
end


%% Set parameters and define variables

Xorig = X;  % save original lesion matrix
res_all.Xorig = Xorig;

Ns = size(X, 1);   % number of subjects


nPC_thresh  = options.nPC_thresh;    % [95% or 99%] percent variance threshold for retaining PCA
l = options.l;

%% Set up variables and do PCA

% exclude voxels with overlapping less than thr_sbj_overlap
idx_voxels  = find(nansum(X, 1) >= options.thr_sbj_overlap);
X           = X(:, idx_voxels);

res_all.idx_voxels = idx_voxels;

% exclude subjects with lesion outsize the threshold
lesion_size         = sum(X, 2);

idx_sbj = ((lesion_size >= options.thr_les_size(1)) & ...
    (lesion_size <= options.thr_les_size(2)));

X = X(idx_sbj, :);
Y = Y(idx_sbj);

res_all.lesion_size = lesion_size;
res_all.idx_sbj     = idx_sbj;

% exclude patients with missing data for DV
ind_y   = find(~isnan(Y));
Y       = Y(ind_y);
X       = double(X(ind_y, :));

% Remove NaN values in X
X(isnan(X) == 1)    = 0; % remove NANs from X

res_all.X = X;

%% Run PCA on the different damage measures

[coeff, score, ~, ~, latent] = pca(X); % lesion masks
res_all.coeff = coeff;

% Get percent variance explained by each component
pve = cumsum(latent);    % variance explained by lesion components

% Find number of PCs necessary to satisfy pve_thresh
nPC = find(pve >= nPC_thresh, 1); % number of components retained
res_all.nPC = nPC;

%% Run ridge regressions with leave-one-out lambda optimization and cross-validation

Xpc = score(:, 1 : nPC);        % lesion predictors
res_all.Xpc = Xpc;

% Normalization
mX = mean(Xpc(:));  % mean across the whole matrix
sX = std(Xpc(:));   % standard deviation of the whole matrix
Xpc = (Xpc - mX) ./ sX;

res_all.Xpc_zscore  = Xpc;
res_all.mX          = mX;
res_all.sX          = sX;

coeff_new = coeff(:, 1 : nPC) ./ sX;
res_all.coeff_new = coeff_new;

if options.zscoreY
    Y   = zscore(Y, [], 1);
end
res_all.mY = mean(Y);
res_all.sY = std(Y);

res_all.Y = Y;

Ypred = nan(length(Y), length(l));

beta_all = nan(nPC, length(Y), length(l));

RSS   = zeros(length(Y), length(l));
tic
for il = 1 : length(l)
% parfor il = 1 : length(l)

    lambda = l(il);
    Yp = nan(length(Y), 1);
    beta_il = nan(nPC, length(Y));
    for is = 1 : length(Y)
    
        idx_train = setdiff(1 : length(Y), is);
        X_train   = Xpc(idx_train, :);
        Y_train   = Y(idx_train);

        X_test    = Xpc(is, :);
        
        betas = pinv(X_train' * X_train + lambda * ...
            eye(nPC)) * X_train' * Y_train;

        Yp(is) = X_test * betas;
        
        beta_il(:, is) = betas;
        
        
    end
    beta_all(:, :, il) = beta_il;
    Ypred(:, il) = Yp;
    
    RSS(:, il) = abs(Ypred(:, il) - Y);
end

RMSD = sqrt(mean(RSS .^ 2, 1));

[RMSD_opt, l_ind] = min(RMSD);
lambda_opt        = l(l_ind);
disp(['lambda_opt = ', num2str(lambda_opt)])

figure;
plot(log10(l), RMSD, '-k')
hold on
plot(log10(lambda_opt), RMSD_opt, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k')
hold off
xlabel('lambda')
ylabel('RMSD')

toc

beta_opt = squeeze(beta_all(:, :, l_ind));
res_all.beta_ave = nanmean(beta_opt, 2);

RR2 = 1 - sum((Y - Ypred(:, l_ind)) .^2) / sum((Y - mean(Y)).^2);
res_all.RR2 = RR2;
res_all.Ypred = Ypred(:, l_ind);

disp('Model estimate: done!')
disp('beta_ave = ')
disp(res_all.beta_ave)
disp('------')
disp(['R2 = ', num2str(RR2)])

% Model significance test

Nperm       = options.Nperm;
RR2_perm    = zeros(Nperm, 1);
beta_perm   = nan(nPC, Nperm);

res_all.Nperm = Nperm;

disp('-------------------------------------------------------------------')
disp('Permutation starts...')
disp('R2 perm = ')

parfor ip = 1 : Nperm
    beta_all        = nan(nPC, length(Y));
    ind_perm        = randperm(length(Y));
    Y_perm          = Y(ind_perm);
    Ypred_perm      = nan(length(Y), 1);
    
    RSS_perm        = zeros(length(Y), length(l));
    
    for il = 1 : length(l)
        lambda = l(il);
        Yp = nan(length(Y), 1);
        beta_il = nan(nPC, length(Y));
        for is = 1 : length(Y)
    
            idx_train = setdiff(1 : length(Y), is);
            X_train   = Xpc(idx_train, :);
            Y_train   = Y_perm(idx_train);

            X_test    = Xpc(is, :);
        
            betas = pinv(X_train' * X_train + lambda * ...
                eye(nPC)) * X_train' * Y_train;

            Yp(is) = X_test * betas;
        
            beta_il(:, is) = betas;
        end
        
        beta_all(:, :, il) = beta_il;
        Ypred_perm(:, il)  = Yp;
    
        RSS_perm(:, il) = abs(Ypred_perm(:, il) - Y_perm);
    end
    
    RMSD = sqrt(mean(RSS_perm .^ 2, 1));

    [~, l_ind] = min(RMSD);
    
    beta_opt = squeeze(beta_all(:, :, l_ind));
    beta_perm(:, ip) = nanmean(beta_opt, 2);

    RR2_perm(ip) = 1 - sum((Y_perm - Ypred_perm(:, l_ind)) .^2) / ...
        sum((Y_perm - mean(Y_perm)).^2);
    disp(num2str(RR2_perm(ip)))
    
end

res_all.beta_perm   = beta_perm;
res_all.RR2_perm    = RR2_perm;

