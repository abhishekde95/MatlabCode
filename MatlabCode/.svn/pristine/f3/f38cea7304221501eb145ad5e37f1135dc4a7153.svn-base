function [nextX, datastruct] = computeNextStim_ALalgorithm_datastruct(datastruct)

tic;

% update data structure for the rest
datastruct.nstim = length(datastruct.r);
thTrial = length(datastruct.r);

if ((rem(thTrial, 10)==0) &&(datastruct.detH>datastruct.thrsh_detH))
    % optimize hyperparameters with analytic form
    ovrscl_1 = mean(datastruct.r)/2; % overall scale
    lngthscl_1 = max(max(datastruct.support))-min(min(datastruct.support))/2; % variance
    prs0 = [mean(datastruct.r); ovrscl_1; lngthscl_1];
    datastruct.K = abs(prs0(2))*exp(-.5/abs(prs0(3)).*datastruct.norm_mat);
    datastruct.muf = abs(prs0(1));
    
    [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
    datastruct.detH=detH;
    datastruct.prs = prs;
    datastruct.finit = fmapFinal;
    datastruct.muf = prs(1);

else
    prs = datastruct.prs;
    datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
    datastruct.muf = abs(prs(1));
    [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct);
    
    datastruct.finit = fmapFinal;
end

datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
[predictiveMean, predictiveVar, nextX] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, datastruct.fvar_logexp1);

toc;
%%
figure(300); clf;
subplot(223);
qz = datastruct.fmean_logexp1(predictiveMean, sqrt(predictiveVar));
npts = 25;
support_x = 1:npts;
support_y = 1:npts;
[xx, yy] = meshgrid(support_x, support_y);

contour(xx, yy, reshape(qz, length(xx), []), 5); axis image; axis xy; title('log-exp-g');
hold on;
plot(datastruct.x(1:end, 1), datastruct.x(1:end,2), 'mo', 'LineWidth',2,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor',[.2 0.2 .6],...
    'MarkerSize',8);
subplot(224);


datastruct.mse(thTrial) = norm(datastruct.lambdatrue - qz);
plot(datastruct.numInitData:thTrial, datastruct.mse(datastruct.numInitData:thTrial), 'ro');

pause(0.01);