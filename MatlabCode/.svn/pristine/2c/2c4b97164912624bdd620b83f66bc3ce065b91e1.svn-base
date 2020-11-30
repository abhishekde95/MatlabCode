function [wt_est] = weights_test(stimuli, resp_act, preproc)
%% Test maximizing the weights 

%% Preproceess
if preproc
    stimuli = zscore(stimuli,0,1);
    resp_act = resp_act - mean(resp_act);
end

%% Define stimuli
n = size(stimuli, 2);
nstim = size(stimuli, 2);
%% Initial Guess
% wts_0 = [15 5 1 1 1, 2, 3];
wts_0 = ones(1,n)*2;
% wts_0 = [10,10,3,1,1,1,2,-1,1,1];

%% Numerically minimize MSE + LASSO
lambda = 0; % Set to 0 to remove LASSO regularization
func_ = @(weights)STS_mu(weights, stimuli, lambda, resp_act); % Function 
% tic
% fminunc method
% options = optimoptions('fminunc', 'Display', 'final', ...
%     'MaxIterations', 10e10, 'Algorithm', 'quasi-newton', ...
%     'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 10e10);
% wt_est = fminunc(func_, wts_0, options);
% toc
% fminsearch method
% tic
options = optimset('Display', 'final', 'MaxIter', 10e10,...
    'MaxFunEvals', 10e10);
wt_est = fminsearch(func_, wts_0, options);
% toc
% Global search method
% gs = GlobalSearch('Display', 'iter');
% opts = optimoptions(@fmincon,'MaxIterations', 1e10, ...
%     'MaxFunctionEvaluations', 1e5);
% prb = createOptimProblem('fmincon', 'x0', wts_0, 'objective', func_, ...
%     'options', opts);
% wt_est = run(gs, prb);

%% Calculate STA with optimized weights
wt_opt_norm = wt_est./max(wt_est);
resp_est = wt_est*stimuli';
STS = resp_est*stimuli;
STA = STS./nstim;
wt_est;
STA;

%% Cost function
function [mu] = STS_mu(weights, stimuli, lambda, resp_act)
% MSE + LASSO
mu = sum((resp_act-(weights*stimuli')).^2, 'all')/length(resp_act) + lambda*norm(weights, 1);
end
end
