function  [neglogev, prs, hmapFinal, aFinal, WFinal, sqrtLFinal, H] = setHyperparam_usingAnalyticForm_nse(prs0, datastruct)

% set bounds on estimated parameters
% rMax = max(datastruct.r);
% xMax = max(max(datastruct.x) - min(datastruct.x));
% 
% muf = [-1, 2*log(rMax)];
% alpha = [1e-3, 2*log(rMax)]; 
% gamma = [1e-6, 2*xMax]; 
% muf = [-10, 10];
% alpha = [1e-3, 10]; 
% gamma = [1e-3, 10]; 

muf = [-10, 10];
alpha = [-10, 10]; % log alpha
gamma = [-10, 10]; % log gamma
nsevar = [-10, 10]; % log nsevar

% LB = [muf(1); alpha(1); gamma(1)];
% UB = [muf(2); alpha(2); gamma(2)]; 

LB = [muf(1); alpha(1); gamma(1); nsevar(1)];
UB = [muf(2); alpha(2); gamma(2); nsevar(2)]; 

lfun= @(p)computeMarginal_usingAnalyticForm_nse(p, datastruct);

opts = optimset('Display','iter','TolFun',1e-10, 'algorithm', 'active-set', 'maxIter',1e10, 'MaxFunEvals', 5*1e3, 'TolX', 1e-10 );
% ------ Optimize evidence --------------------------------------
[prs, fval, exitflag, output, lambda, grad, hessian] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

[neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal, H] = lfun(prs);
%  [neglogev, hmap, a, W, sqrtL, H]
% fprintf('E.A. is done');

% HessCheck(lfun,prs0);
% HessCheck_Elts(lfun,[1 1],2*prs0);
% DerivCheck(lfun,prs0);
% DerivCheck_Elts(lfun,1,prs0);