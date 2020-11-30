function [h,p] = tmp_calcquadSSE(final_model,RHO,THETA,not_oog_idx)
% calculates sum of squared errors(product of likelihoods), sum of likelihoods, AIC and BIC
RHOtmp = RHO(not_oog_idx);
THETAtmp = THETA(not_oog_idx);
[~,idx] = sort(THETAtmp); % sorting THETA in ascending order
RHOtmp = RHOtmp(idx);
THETAtmp = THETAtmp(idx);
A = final_model(1); B = final_model(2); C = final_model(3); D = final_model(4); E = final_model(5); 
pred = [];
p = [A*cos(THETAtmp).^2+B*sin(THETAtmp).^2+C*(cos(THETAtmp).*sin(THETAtmp)) D*cos(THETAtmp)+E*sin(THETAtmp) -1*ones(numel(THETAtmp),1)];
for kk = 1:size(p,1)
    r = max(roots(p(kk,:)));
    pred = [pred; r];
end
LOOGtmp= pred<0;
resid = log(pred) - log(RHOtmp);

% fvalRob = tmp_calcerror_wo_oog(pred',RHOtmp,2);
% Wald-Wolfowitz Runs test on residuals
if any(LOOGtmp)
    resid(LOOGtmp) = 1;
end
[h,p] = runstest(resid(~LOOGtmp),0,'Alpha',0.05);
end
