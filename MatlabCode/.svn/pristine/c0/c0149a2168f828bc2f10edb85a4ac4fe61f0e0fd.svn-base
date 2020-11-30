%% to test variance-based active learning algorithm
%% define input support and true function

clear;
clc;
% clf;

maxrep = 1;
HowmanyData = 500;
mse = zeros(HowmanyData, maxrep);
% load fvar_logexp1_cubic.mat;
% load fmean_logexp1_cubic.mat;
load fvar_logexp1_lin.mat;
load fmean_logexp1_lin.mat;

fvar_logexp1 = fvar_logexp1_lin;
fmean_logexp1 = fmean_logexp1_lin;

% true lambda
npts = 25;
support_x = 1:npts;
support_y = 1:npts;
[xx, yy] = meshgrid(support_x, support_y);

support = [xx(:) yy(:)];
datastruct.norm_mat_support = form_normMat(support, support);  % squared distance

muvec1 = [npts*0.05 npts*0.05];
muvec2 = [npts*0.5 npts*0.9];
muvec3 = [npts*0.8 npts*0.4];

invC1 = diag([1/npts 1/npts]);
invC2 = diag([4/npts 4/npts]);
invC3 = diag([4/npts 4/npts]);

intensity = @(t) 90*diag(exp(-0.5*(bsxfun(@minus, t, muvec1))*invC1*(bsxfun(@minus, t, muvec1))')) + 80*diag(exp(-0.5*(bsxfun(@minus, t, muvec2))*invC2*(bsxfun(@minus, t, muvec2))')) ...
    + 60*diag(exp(-0.5*(bsxfun(@minus, t, muvec3))*invC3*(bsxfun(@minus, t, muvec3))')) ...
    + sum(log(exp(sin(t*.2.*pi./50).*4-1) + 1), 2);

lambda_true =intensity(support);

g = @(t) log(exp(t)+1);
ginv = @(t) log(exp(t)-1);

figure(300);
subplot(221); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; title('true firing map');
subplot(222); plot(lambda_true); title('true firing map in 1D');

%%
for nrep = 1:maxrep
    
    nrep
    
    %% generate initial stimuli from latin hypercube design
    
    numb_cn = 2;
    numb_tr = numb_cn*10;
%     numb_tr = 500;
    
    % initial points from LH design
    stim_range_min = min(support);
    stim_range_max = max(support);
    x_lh = lhsdesign(numb_tr, numb_cn);
    x = bsxfun(@plus, bsxfun(@times, (stim_range_max - stim_range_min), x_lh), stim_range_min);
    
    lam_x = g(intensity(x));
    r = poissrnd(lam_x);
    r(r==0) = 0.5;
    
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
    datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
    param = abs(prs0(2:end));
    datastruct.K = param(1)*exp(-.5/param(2).*datastruct.norm_mat);% covariance matrix
    
    [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0] = computeFmapAndUpdateTheta(prs0, datastruct);
    
    datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
    datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
    
    %% make prediction
    
    %     [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
    
    %%  Active learning
    
    datastruct.finit = fmapFinal;
    count=1;
    thrsh_detH = 0.005;
    detH = 1;
    
    while(count<=HowmanyData)
        [nrep count]
        
        tic;
        %  predict f given data points
        [predictiveMean, predictiveVar, xNext, idxNext, predictiveCov] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, fvar_logexp1);
        
        gfmean = fmean_logexp1(predictiveMean, sqrt(predictiveVar));
        mse(count, nrep) = norm(lambda_true - gfmean);
        
        % do experiment
        new_r = poissrnd(g(intensity(xNext)));
        
        datastruct.x = [datastruct.x; xNext];
        datastruct.r = [datastruct.r; new_r];
        datastruct.finit = [datastruct.finit; ginv(new_r+eps)];
        datastruct.nstim = (datastruct.nstim+1);
        
        % whenever update x, update norm_mat
        sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
        datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
        %         datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);
        
        if (count~=1) && (rem(count, 10)==0) &&(detH>thrsh_detH)
            % optimize hyperparameters with analytic form
            ovrscl_1 = mean(datastruct.r)/2; % overall scale
            lngthscl_1 = max(max(support))-min(min(support))/2; % variance
            prs0 = [mean(datastruct.r); ovrscl_1; lngthscl_1];
            datastruct.K = abs(prs0(2))*exp(-.5/abs(prs0(3)).*datastruct.norm_mat);
            datastruct.muf = abs(prs0(1));
            [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
        else
            datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
            datastruct.muf = abs(prs(1));
            [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct);
        end
        toc;
        
        %%
        
        datastruct.finit = fmapFinal;
        normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
        datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
        datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
        
        %%
        figure(300); clf;
        subplot(221); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; title('true firing map');
        subplot(223);
        [predictiveMean_onGrid, predictiveVar_onGrid, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, fvar_logexp1);
        
         qz = fmean_logexp1(predictiveMean_onGrid, sqrt(predictiveVar_onGrid));
        
        contour(xx, yy, reshape(qz, length(xx), []), 5); title('log-exp-g');
        hold on;
        plot(datastruct.x(1:end, 1), datastruct.x(1:end,2), 'mo', 'LineWidth',2,...
            'MarkerEdgeColor','y',...
            'MarkerFaceColor',[.2 0.2 .6],...
            'MarkerSize',8);
        subplot(224);
        plot(1:count, mse(1:count, nrep), 'ro');
        
        datastruct.support = support;
        count = count+1;
        disp('br-logexp');
        %         disp('mi-logexp');
        pause(0.5);
    end
    
end

