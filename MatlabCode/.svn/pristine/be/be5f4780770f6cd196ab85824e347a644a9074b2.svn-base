function  [prs, hmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta_nse(prs0, datastruct)

threshEvDiff = 1e-1;

count = 1;
maxIter = 10;

% 1. given theta0, form K
[neglogev0, hmapR0, a0, W0, sqrtL0]  = computeFmap_nse(datastruct);
datastruct.L = sqrtL0.^2;
datastruct.Lminit = datastruct.L*hmapR0 + a0;

datastruct.neglogev = neglogev0;

while(count<maxIter)
    
    % 0. from theta0, optimize hyperparameters using analytic form
    [neglogev1, prs, hmapFinal, aFinal, WFinal, sqrtLFinal, Hfinal, detH] = updateThetaGivenL_nse(prs0, datastruct);
    datastruct.hinit =  hmapFinal;
    datastruct.H = Hfinal; % H = K + nsevar I
    datastruct.mu = abs(prs(1)); 
    
    % 1. given theta0, form H
    [neglogev0, hmapR0, a0, W0, sqrtL0]  = computeFmap_nse(datastruct);
    
    evDiff = datastruct.neglogev - neglogev0; % should be positive, if it goes to the right direction
    
    if (abs(evDiff))<= threshEvDiff
        %         neglog = datastruct.neglog0;
        return;
    elseif evDiff>0
        prs0 = prs;
        datastruct.neglogev = neglogev0;
        datastruct.L = sqrtL0.^2;
        datastruct.Lminit = datastruct.L*hmapR0 + a0;
    else
        prs0 = (prs0+prs)./2;
%         datastruct.finit = fmapR0;
    end
    
    count = count +1; 
    
end