function  [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct)

threshEvDiff = 0.1;

count = 1;
maxIter = 10;

% 1. given theta0, form K
[neglogev0, fmapR0, a0, W0, sqrtL0]  = updateFmapGivenK(datastruct);
datastruct.L = sqrtL0.^2;
datastruct.Lminit = datastruct.L*fmapR0 + a0;
datastruct.neglogev = neglogev0;
% datastruct.finit =  fmapR0;

while(count<maxIter)
    
    % 0. from theta0, optimize theta using analytic form
    [neglogev1, prs, fmapFinal, aFinal, WFinal, sqrtLFinal, Kfinal, detH] = updateThetaGivenL(prs0, datastruct);
    datastruct.finit =  fmapFinal;
    datastruct.K = Kfinal;
    datastruct.mu = abs(prs(1));
    
    % 1. given theta0, form K
    [neglogev0, fmapR0, a0, W0, sqrtL0]  = updateFmapGivenK(datastruct);
    
    evDiff = datastruct.neglogev - neglogev0; % should be positive, if it goes to the right direction
    
    if (abs(evDiff))<= threshEvDiff
        %         neglog = datastruct.neglog0;
        return;
    elseif evDiff>0
        prs0 = prs;
        datastruct.neglogev = neglogev0;
        datastruct.L = sqrtL0.^2;
        datastruct.Lminit = datastruct.L*fmapR0 + a0;
    else
        prs0 = (prs0+prs)./2;
%         datastruct.finit = fmapR0;
    end
    
    count = count +1;
    
end
