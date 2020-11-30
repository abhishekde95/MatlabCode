function [f, fmap, a, W, sqrtL, K]  = updateHyperparam_wAnalyticForm(p, datastruct)

% 1. given theta, form K
muf = abs(p(1));
param = abs(p(2:end));
norm_mat = datastruct.norm_mat; 
K = param(1)*exp(-.5/param(2).*norm_mat);

% 2. given K, find fmap
Lm_init = datastruct.Lminit;
sqrtL = sqrt(datastruct.L);

B = eye(size(sqrtL)) + sqrtL*K*sqrtL;

W = chol(B, 'lower');

M = W\sqrtL;
% sqrtLinvBsqrtL = sqrtL*(W'\(W\sqrtL));
sqrtLinvBsqrtL = M'*M;

a = Lm_init - sqrtLinvBsqrtL*(K*Lm_init + muf);

fmap = K*a + muf;

% 3. given fmap, compute neglogev
f = computeLogevidence(fmap, W, muf, a, datastruct);