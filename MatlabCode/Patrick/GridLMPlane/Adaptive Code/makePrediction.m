function [predictiveMean, predictiveVar, xNext, idxNext, predictiveCov] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, fvar_logexp1, whichMethod)
% [predictiveMean, predictiveVar, nextX] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, datastruct.fvar_logexp1, datastruct.whichMethod);


param = prs(2:end);
% Kstar = formK(datastruct.support, datastruct.x, param);
% Kstarstar = formK(datastruct.support, datastruct.support, param);
Kstar = datastruct.Kstar;
norm_mat_support = datastruct.norm_mat_support;
Kstarstar = param(1)*exp(-.5/param(2).*norm_mat_support);

muf = prs(1);
predictiveMean = muf + Kstar*aFinal;
predictiveCov = (Kstarstar - Kstar*sqrtLFinal*(WFinal'\(WFinal\sqrtLFinal))*Kstar');

predictiveVar = diag(predictiveCov);

%% uncertainty sampling:

% [g, dg] = logexp1(predictiveMean);
% var_lambda = dg.*sqrt(predictiveVar);
%
% idx = find(var_lambda == max(var_lambda));
% idxNext = idx(floor(rand*length(idx))+1);
% xNext = datastruct.support(idxNext,:);

%% approx. Bayes risk

if whichMethod ==1
    
    [g, dg] = logexp1(predictiveMean);
    sleng = size(predictiveCov,1);
    mu_t1 = predictiveMean;
    % sig_t1 = zeros(sleng, sleng);
    %
    % for j=1:sleng
    %     sig_t1(:,j) = predictiveVar - (predictiveCov(:,j).^2.*dg(j).^2./g(j))./(1+predictiveVar(j).*dg(j).^2./g(j));
    % end
    
    dGtimesG = dg.^2./g;
    CovTimesDgbyG = bsxfun(@times, predictiveCov.^2, dGtimesG');
    
    VardGtimesG = 1+ predictiveVar.*dGtimesG;
    secondTrm = bsxfun(@times, CovTimesDgbyG, 1./VardGtimesG');
    sig_t1 = bsxfun(@minus, predictiveVar, secondTrm);
       
%     keyboard;
   
    var_lambda = fvar_logexp1(repmat(mu_t1, 1, sleng), sqrt(sig_t1));
    
    sumvarlam = sum(var_lambda);
    % idx = find(sumvarlam == min(sumvarlam));
    idx = find(sumvarlam == min(sumvarlam));
    
    if isempty(idx)
        idxNext = floor(rand*length(datastruct.support))+1;
    else
        idxNext = idx(floor(rand*length(idx))+1);
    end
    xNext = datastruct.support(idxNext,:);
    
%% mutual information    
else % whichMethod ==2
    
%     load datastruct2Dtable;
%     
%     fxx = datastruct2Dtable.fxx;
%     sigxx = datastruct2Dtable.sigxx;
%     val = datastruct2Dtable.val;
    
    [xi, yi] = meshgrid(predictiveMean, sqrt(predictiveVar));
    mutualInfo = interp2(datastruct.fxx, datastruct.sigxx, datastruct.val, xi, yi, '*linear');
    MI = diag(mutualInfo);
    
    idx = find(MI == max(MI));
    if size(idx, 1) == 0
        % in case predictiveMean is all the same constant
        idxNext = randi(length(datastruct.support), 1, 1);
    else
        idxNext = idx(floor(rand*length(idx))+1);
    end
    
    xNext = datastruct.support(idxNext,:);
    
end
