
%% to test our new algorithm (26th of March, 2012)
% 2D example! is is really better in high-D?

%% define input support and true function

clear;
clc;
% clf;

maxrep = 1;
% mse_br = zeros(maxrep, 1);
HowmanyData = 500;
mse = zeros(HowmanyData, maxrep);
postent = zeros(HowmanyData, maxrep);
totVar = zeros(HowmanyData, maxrep);
load fvar_logexp1_cubic.mat;
load fmean_logexp1_cubic.mat;
load gfmean_2Dalldata.mat;

%%
for nrep = 1:maxrep
    
    nrep
    %%
    [num, txt, raw] = xlsread('statmat2.xls');
    
    xraw = num(:,1:2);
    spkraw = num(:,3);
    
    %%
    
    [support, i_uniq] = unique(xraw, 'rows');
    spk = spkraw(i_uniq, :);
    
    [sortX,i] = sortrows(xraw);
%     support = sortX(1:3:end, :);
    spksort = reshape(spkraw(i),3,[])';
    
    %% generate initial stimuli from latin hypercube design
    
    numb_cn = 2;
    numb_tr = numb_cn*10;
    
    % the max and min stimulus range
%     stim_range_min = min(support);
%     stim_range_max = max(support);
%     x_lh = lhsdesign(numb_tr, numb_cn);
%     x = bsxfun(@plus, bsxfun(@times, (stim_range_max - stim_range_min), x_lh), stim_range_min);
     
%     idx_stim = lhsdesign(numb_tr, numb_cn);
%     idx_x = bsxfun(@times, length(support), idx_stim);
%     idx_x = round(idx_x);
    
%     zi = interp2(support(:,1), support(:,2), spksort(:,1), x(:,1), x(:,2));
    idx_x = round(length(support)*rand(numb_tr, 1));
    x = support(idx_x, :);
    
    r = spk(idx_x, :);
    r(r==0) = 0.5;

    g = @(t) log(exp(t)+1);
    ginv = @(t) log(exp(t)-1);
    
    %% Given data, optimize theta and find fmap
    
    % data structure
    datastruct.support = support;
    datastruct.r = r;
    datastruct.x = x;
    datastruct.nstim = numb_tr;
    datastruct.ndim = numb_cn;
    datastruct.finit = ginv(r+0.1);
    datastruct.g = g;
    datastruct.ginv = ginv;
    
    % initial hyperparameters
    datastruct.muf = mean(r);
    ovrscl_1 = mean(r)/2; % overall scale
    lngthscl_1 = max(max(support))-min(min(support))/2; % variance
    prs0 = [datastruct.muf; ovrscl_1; lngthscl_1];
    
    % optimization with analytic form
    datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);   
%     param = abs(prs0(2:end));
%     datastruct.K = param(1)*exp(-.5/param(2).*datastruct.norm_mat);
    
    [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0] = updateFmapHyperparam_main(prs0, datastruct);
    
    %% make prediction
    
%     [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
    
    %%  Active learning
    
    datastruct.finit = fmapFinal;
    count=1;
    
    while(count<=HowmanyData)
        [nrep count]
        
        tic;
        %  predict f given data points
        [predictiveMean, predictiveVar, xNext, idxNext, predictiveCov] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
%         totVar(count, nrep) = sum(fvar_logexp1_cubic(predictiveMean, sqrt(predictiveVar)));
        
        gfmean = fmean_logexp1_cubic(predictiveMean, sqrt(predictiveVar));
        mse(count, nrep) = norm(gfmean_2Dalldata - gfmean);
%         postent(count, nrep) = logdetns(predictiveCov);

        % do experiment
        new_r = spk(idxNext, :);
        
        datastruct.x = [datastruct.x; xNext];
        datastruct.r = [datastruct.r; new_r];
        datastruct.finit = [datastruct.finit; ginv(new_r+eps)];
        datastruct.nstim = (datastruct.nstim+1);
        
        % whenever update x, update norm_mat
        datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);
%         datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
        
        if (count~=1) && (rem(count, 10)==0)
            % optimize hyperparameters with analytic form
            ovrscl_1 = mean(datastruct.r)/2; % overall scale
            lngthscl_1 = max(max(support))-min(min(support))/2; % variance
            prs0 = [mean(datastruct.r); ovrscl_1; lngthscl_1];
            [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0] = updateFmapHyperparam_main(prs0, datastruct);
        else
            [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmap(prs, datastruct);
        end
        toc;
        
        datastruct.finit = fmapFinal;
               
        
        figure(300); clf;
        subplot(221); 
        lengAxis = 20;
        xaxis = linspace(min(support(:,1)), max(support(:,1)), lengAxis);
        yaxis = linspace(min(support(:,2)), max(support(:,2)), lengAxis);
        
        [xx, yy] = meshgrid(xaxis, yaxis);
        datastruct.support = [xx(:) yy(:)];        
        [predictiveMean_onGrid, predictiveVar_onGrid, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
        qz = fmean_logexp1_cubic(predictiveMean_onGrid, sqrt(predictiveVar_onGrid));
        
        contour(xx, yy, reshape(qz, length(xx), []), 5); title('log-exp-g');
        hold on;
        plot(datastruct.x(numb_tr+1:end, 1), datastruct.x(numb_tr+1:end,2), 'mo', 'LineWidth',2,...
            'MarkerEdgeColor','y',...
            'MarkerFaceColor',[.2 0.2 .6],...
            'MarkerSize',8);
        subplot(222);
        plot(1:count, mse(1:count, nrep), 'ro');
        
        datastruct.support = support;
        count = count+1;
        disp('br-logexp');
%         disp('mi-logexp');

    end
    
end

