function [h,p,nInF,fvalRob,fvalOLS,Rlin,adjR2] = calclinSSE(final_model,RHO,THETA,not_oog_idx)
% calculates sum of squared errors(product of likelihoods), sum of
% likelihoods, AIC and BIC, and runs test for the residuals

RHOtmp = RHO(not_oog_idx);
THETAtmp = THETA(not_oog_idx);
[~,idx2] = sort(THETAtmp); % sorting THETA in ascending order
RHOtmp = RHOtmp(idx2);
THETAtmp = THETAtmp(idx2);
pred = 1./(final_model*[cos(THETAtmp'); sin(THETAtmp')]);
LOOGtmp= pred<0;
resid = log(pred) - log(RHOtmp');
nInF = sum(LOOGtmp);
Rlin = corr(log(pred'),log(RHOtmp));
fvalOLS = calcSSE(pred,RHOtmp,1);
fvalRob = calcSSE(pred,RHOtmp,2); % Equivalent to SSE
fvalRob_Total = calcSSE(mean(RHOtmp),RHOtmp,2); % Equivalent to SST
SST = fvalRob_Total;
SSE = fvalRob;
R2 = (SST-SSE)/(SST-SSE+1.2076*SSE);
n = numel(not_oog_idx);
adjR2 = 1-((1-R2)*(n-1)/(n-2));

% Wald-Wolfowitz Runs test on residuals
if any(LOOGtmp)
    resid(LOOGtmp) = 1;
end
[h,p] = runstest(resid,0,'Alpha',0.05);

end

