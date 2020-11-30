function  [neglogev, prs, fmapFinal, aFinal, WFinal, sqrtLFinal, Kfinal, detH] = updateThetaGivenL(prs0, datastruct)

% set bounds on estimated parameters
muf = [1e-3, 1e3];
alpha = [1e-3, 1e3];
gamma = [1e-3, 1e3];

LB = [muf(1); alpha(1); gamma(1)];
UB = [muf(2); alpha(2); gamma(2)]; 

lfun= @(p)updateHyperparam_wAnalyticForm(p, datastruct);

% opts = optimset('Display','iter','TolFun',1e-6, 'algorithm', 'active-set', 'gradobj', 'on', 'maxIter', 10, 'maxFunEvals', 1e3);
opts = optimset('Display','iter','TolFun',1e-6, 'algorithm', 'active-set', 'maxIter', 10, 'maxFunEvals', 1e3);

% ------ Optimize evidence --------------------------------------
[prs, fval, exitflag, output, lambda, grad, hessian] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

[neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal, Kfinal] = lfun(prs);
%  [f, df, fmap, a, W, sqrtL, K]  = updateHyperparam_wAnalyticForm(p, datastruct)
detH = det(hessian);

% fprintf('E.A. is done');

% HessCheck(lfun,prs0);
% HessCheck_Elts(lfun,[1 1],2*prs0);
% DerivCheck(lfun,prs0);
% DerivCheck_Elts(lfun,1,prs0);