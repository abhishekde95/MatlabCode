function [mean_lambda, xNext, idxNext, var_lambda] = makePrediction_nse(prs, datastruct, aFinal, WFinal, sqrtLFinal)

muf = prs(1);
param = prs(2:3);
nv = prs(end);

Kstar = datastruct.Kstar;
norm_mat_support = datastruct.norm_mat_support;
Kstarstar = param(1)*exp(-.5/param(2).*norm_mat_support);

predictiveMean = muf + Kstar*aFinal;

M = WFinal\(sqrtLFinal*Kstar');

predictiveCov = Kstarstar - M'*M;

% predictiveCov = (Kstarstar - Kstar*sqrtLFinal*(WFinal'\(WFinal\sqrtLFinal))*Kstar');

predictiveVar = diag(predictiveCov);

mean_lambda =  exp(predictiveMean + 0.5*(predictiveVar+nv));
%% uncertainty

var_lambda = exp(nv).*(exp(predictiveVar)-1).*exp(2*predictiveMean + predictiveVar);

idx = find(var_lambda == max(var_lambda));
idxNext = idx(floor(rand*length(idx))+1);
xNext = datastruct.support(idxNext,:);
