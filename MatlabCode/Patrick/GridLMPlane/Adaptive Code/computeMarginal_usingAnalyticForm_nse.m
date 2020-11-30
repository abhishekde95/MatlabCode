function [neglogev, hmap, a, W, sqrtL, H]  = computeMarginal_usingAnalyticForm_nse(p, datastruct)

% 1. given theta, form H
muf = abs(p(1));
param = abs(p(2:3));
nsevar = abs(p(end));

norm_mat = datastruct.norm_mat; 
K = param(1)*exp(-.5/param(2).*norm_mat);
nsevarI = nsevar*eye(size(norm_mat)); 
H = K + nsevarI;

% 2. given K, find fmap
Lm_init = datastruct.Lminit;
sqrtL = sqrt(datastruct.L);

B = eye(size(sqrtL)) + sqrtL*H*sqrtL;

W = chol(B, 'lower');

M = W\sqrtL;
sqrtLinvBsqrtL = M'*M;
a = Lm_init - sqrtLinvBsqrtL*(H*Lm_init + muf);

hmap = H*a + muf;

% 3. given fmap, compute neglogev
neglogev = computeLogevidence_nse(hmap, W, muf, a, datastruct.r);