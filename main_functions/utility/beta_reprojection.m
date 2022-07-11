% Function to reproject the PC-beta to the original voxel
% X = Xpc * coeff' + repmat(mean(X, 1), [Ns 1])
% => Xpc = (X - repmat(mean(X, 1), [Ns 1])) * coeff
% zXpc = (Xpc - repmat(mean(Xpc(:)), [Ns nPC])) ./ std(Xpc(:))
%      = ((X - repmat(mean(X, 1), [Ns 1])) * coeff - repmat(mean(Xpc(:)), [Ns nPC])) ./ std(Xpc(:))
% => xXpx * beta = (X*coeff*beta)./std(Xpc) + 
%                  - (repmat(mean(Xpc(:)),[Ns nPC])*coeff*beta)./std(Xpc(:))
%                  - (repmat(mean(Xpc(:)), [Ns nPC]) * beta) ./ std(Xpc(:))
% => beta_voxel = (coeff*beta)./std(Xpc) + 
%                  - X^-1 * (repmat(mean(Xpc(:)),[Ns nPC])*coeff*beta)./std(Xpc(:))
%                  - X^-1 * (repmat(mean(Xpc(:)), [Ns nPC]) * beta) ./ std(Xpc(:))

function beta_voxel = beta_reprojection(beta, X, coeff, mPC, sPC)

[Ns, n] = size(X);  % Ns = # subject
                    % n = # voxels

nPC = size(coeff, 2);  % # PCs

meanX = repmat(mean(X, 1), [Ns, 1]);
mPC   = repmat(mPC, [Ns, nPC]);     % mean of the PCs-scores

X_inv = pinv(X);

f1 = coeff * beta;
% f2 = X_inv * (meanX * coeff * beta);
f3 = X_inv * (mPC * beta);

% beta_voxel = (f1 - f2 - f3) ./ sPC;
beta_voxel = (f1 - f3) ./ sPC;

% check
if (size(beta_voxel, 1) ~= n) || (size(beta_voxel, 2) ~= 1)
    disp('beta size are wrong!!')
end
