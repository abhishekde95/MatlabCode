function [minmse, mseStep, LARSmse, cvidx] = WNLARS_kfoldcv(X, y, kfold, predictSet, verbose)
%% 
% This function performs kfold cross validation (cv) on LARS outputs to 
% estimate the optimal step size(i.e. feature selection).
% 
% Inputs:
%   'X' [nstim, ncovariates] - 2D array of stimuli data
%   'y' [nstim] - Vector of response
%   'kfold' [int] - Number of folds for k-fold cv
%   'predictSet' [bool] - Exclude a kfold set from LARS and cv
%   'verbose' [bool] - Display which fold is being processed
%
% Outputs:
%   'minmse' [int] - Number of steps corresponding to the minimum MSE
%       (After MSE is averaged across folds)
%   'mseStep' [ncovariates] - Vector of MSE results averaged across cv 
%       folds at each step returned by LARS
%   'LARSmse' [nfolds, ncovariates] - MSE results from each fold of cv for
%       each step of LARS
%   'cvidx' [nstim] - cv kfold labels for each row
%%

% Run cv on all kfold sets if not specified
if nargin < 4; predictSet = False; end
% Do not display kfold progress if not specified
if nargin < 5; verbose = False; end

% Define set of kfold indices to use in cross validation.
if predictSet
    kfoldSet = 1:(kfold-1);
else
    kfoldSet = 1:kfold;
end
    
% Partition data for cross validation
cvidx = crossvalind('Kfold', size(X,1), kfold);

% K-fold cross validation
for ii = kfoldSet
    
    % Define set of indices for LARS for each kfold
    kfoldidx = kfoldSet(kfoldSet~=ii);
    
    if verbose; disp(['Computing k-fold ' num2str(ii)]); end
    
    % Run LARS on training data
    [beta, ~] = WN_LARS(X(ismember(cvidx, kfoldidx),:), ...
        y(ismember(cvidx, kfoldidx)), true, true);
    
    % Ensure test data is centered and normalized appropriately
    Xc = X(cvidx == ii,:) - mean(X(cvidx == ii,:),1);
    Xc = Xc./sqrt(sum(Xc.^2,1));
    yc = y(cvidx == ii)-mean(y(cvidx == ii));
    
    % Compute MSE on testing data
    LARSmse(ii,:) = sum((yc-Xc*beta).^2,1)/size(Xc, 1);
end

% Avg MSE across cv folds
mseStep = mean(LARSmse, 1);

% Find the index of the minimum MSE
[~, minmse] = min(mseStep);