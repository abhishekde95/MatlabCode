function [neglogev, prs, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal, detH] = updateThetaGivenL_nse(prs0, datastruct)

% set bounds on estimated parameters
% these parameters are for f, i.e., log-lambda 
muf = [1e-3, 10]; % mean of f 
alpha = [1e-3, 10]; % overall scale of f
gamma = [1e-3, 1e3]; % smoothness scale of f
nsevar = [1e-6, 10]; % nsevar

LB = [muf(1); alpha(1); gamma(1); nsevar(1)];
UB = [muf(2); alpha(2); gamma(2); nsevar(2)]; 

lfun= @(p)computeMarginal_usingAnalyticForm_nse(p, datastruct);

opts = optimset('Display','iter','TolFun',1e-10, 'algorithm', 'active-set', 'maxIter', 10);

% ------ Optimize evidence --------------------------------------
[prs, fval, exitflag, output, lambda, grad, hessian] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

[neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal] = lfun(prs);
detH = 1/det(hessian); 