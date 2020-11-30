function [nextX, totData] = computeNextStim_ALalgorithm(x, r, stimDim, support,  numInitData, totNumTrials, lambda_true, whichMethod)

tic;
persistent datastruct;

if isempty(datastruct)
    datastruct.x = x;
    datastruct.r = r;
    datastruct.support = support;
    datastruct.ndim = stimDim;
    datastruct.norm_mat_support = form_normMat(support, support);  % squared distance
    
    g = @(t) log(exp(t)+1);
    ginv = @(t) log(exp(t)-1);
    datastruct.g = g;
    datastruct.ginv = ginv;
    datastruct.finit = datastruct.ginv(datastruct.r+0.1);
    
    datastruct.thrsh_detH = 0.01;
    datastruct.detH = 1;
    
    datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
    datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
    
    load fvar_logexp1_lin.mat;
    load fmean_logexp1_lin.mat;
    
    datastruct.fvar_logexp1 = fvar_logexp1_lin;
    datastruct.fmean_logexp1 = fmean_logexp1_lin;
    
    datastruct.whichMethod = whichMethod; 
    
    load datastruct2Dtable;
    
    datastruct.fxx = datastruct2Dtable.fxx;
    datastruct.sigxx = datastruct2Dtable.sigxx;
    datastruct.val = datastruct2Dtable.val;
    
    % delete this for real experiments
    datastruct.mse = zeros(totNumTrials-numInitData,1);
    datastruct.lambdatrue = lambda_true;
    
else
    datastruct.x = x;
    datastruct.r = r;
    sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
    datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
    
    normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
    datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
    datastruct.finit = [datastruct.finit; datastruct.ginv(r(end)+0.1)];
end

% update data structure for the rest
datastruct.nstim = length(datastruct.r);
thTrial = length(datastruct.r);


if ((rem(thTrial, 10)==0) &&(datastruct.detH>datastruct.thrsh_detH))
    % optimize hyperparameters with analytic form
    ovrscl_1 = mean(datastruct.r)/2; % overall scale
    lngthscl_1 = max(max(support))-min(min(support))/2; % variance
    prs0 = [mean(datastruct.r); ovrscl_1; lngthscl_1];
    datastruct.K = abs(prs0(2))*exp(-.5/abs(prs0(3)).*datastruct.norm_mat);
    datastruct.muf = abs(prs0(1));
    neglogev0  = updateFmapGivenK(datastruct);
    
    if thTrial> numInitData
        if datastruct.neglogev<neglogev0
            prs = datastruct.prs;
            datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
            datastruct.muf = abs(prs(1));
            [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct);
            
            datastruct.neglogev = neglogev;
            datastruct.finit = fmapFinal;
        else
            [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
            datastruct.detH=detH;
            datastruct.prs = prs;
            datastruct.finit = fmapFinal;
            datastruct.muf = prs(1);
        end
    else
        [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
        datastruct.detH=detH;
        datastruct.prs = prs;
        datastruct.finit = fmapFinal;
        datastruct.muf = prs(1);
    end
               

else
    prs = datastruct.prs;
    datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
    datastruct.muf = abs(prs(1));
    [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct);
    
    datastruct.neglogev = neglogev; 
    datastruct.finit = fmapFinal;
end

datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
[predictiveMean, predictiveVar, nextX, idxNext, predictiveCov] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, datastruct.fvar_logexp1, datastruct.whichMethod);

toc;

qz = datastruct.fmean_logexp1(predictiveMean, sqrt(predictiveVar));

datastruct.lambFinal = qz; 
datastruct.aFinal = aFinal;
datastruct.WFinal = WFinal;
datastruct.sqrtLFinal = sqrtLFinal;
totData = datastruct;