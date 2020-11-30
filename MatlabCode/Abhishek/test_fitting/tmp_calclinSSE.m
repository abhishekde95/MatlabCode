function [h,p] = tmp_calclinSSE(final_model,RHO,THETA,not_oog_idx)
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

% fvalRob = tmp_calcerror_wo_oog(pred,RHOtmp,2);
% Wald-Wolfowitz Runs test on residuals
if any(LOOGtmp)
    resid(LOOGtmp) = 1;
end
[h,p] = runstest(resid,0,'Alpha',0.05);

end


