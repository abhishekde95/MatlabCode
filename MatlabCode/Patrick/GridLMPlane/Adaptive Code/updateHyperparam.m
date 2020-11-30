function  [neglogev, prs, fmapFinal, aFinal, WFinal, sqrtLFinal, Kfinal] = updateHyperparam(prs0, datastruct)

% set bounds on estimated parameters
muf = [1e-3, 1e3];
alpha = [1e-6, 1e6];
gamma = [1e-6, 1e6];
% l0 =  [1e-3, 5];
% l1 = [1e-6, 1e2];

% muf = [-10, 10];
% alpha = [-10, 10]; % log alpha
% gamma = [-8, 8]; % log gamma
% gamma = [-10, 10]; % log gamma
% gamma = [-10, 10]; % log gamma

% LB = [muf(1); alpha(1); gamma(1)];
% UB = [muf(2); alpha(2); gamma(2)]; 

LB = [muf(1); alpha(1); gamma(1)*ones(length(prs0)-2,1)];
UB = [muf(2); alpha(2); gamma(2)*ones(length(prs0)-2,1)]; 

lfun= @(p)updateHyperparam_wAnalyticForm(p, datastruct);

opts = optimset('Display','iter','TolFun',1e-6, 'algorithm', 'active-set', 'maxIter', 10, 'maxFunEvals', 1e3);
% ------ Optimize evidence --------------------------------------
[prs, fval, exitflag, output, lambda, grad, hessian] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

[neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal, Kfinal] = lfun(prs);

% fprintf('E.A. is done');

% HessCheck(lfun,prs0);
% HessCheck_Elts(lfun,[1 1],2*prs0);
% DerivCheck(lfun,prs0);
% DerivCheck_Elts(lfun,1,prs0);