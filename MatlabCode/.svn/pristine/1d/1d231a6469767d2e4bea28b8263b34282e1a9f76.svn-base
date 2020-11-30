%% to test our new algorithm to Greg's data (10 th of March, 2012)
% new components:
%   1. no inv(K)
%   2. exponential g
%   3. max info criterion
%   4. test purely numerical optimization and mixture of numerical and
%   analytic optimization for hyperparameters.

%% define input support and true function

clear;
clc;

[num, txt, raw] = xlsread('statmat2.xls');

x = num(:,1:2);
spk = num(:,3);

[support_x, i] = unique(x, 'rows');
% support_x = x(i, :);
% spk = spk(i,:);

%%
% [sortX,i] = sortrows(x);
% spksort = reshape(spk(i),3,[])';
% muspk = mean(spksort,2);

%% generate simulated data

numb_cn = 2;
numb_tr = length(x);
r = spk; % to avoid numerical problems

r(r==0) = 0.5; 

%% Given data, optimize theta and find fmap

g = @(t) log(exp(t)+1);
ginv = @(t) log(exp(t)-1);

% data structure
datastruct.support = support_x;
datastruct.r = r;
datastruct.x = x;
datastruct.nstim = numb_tr;
datastruct.ndim = numb_cn;
datastruct.finit = ginv(r);
datastruct.g = g;
datastruct.ginv = ginv;

%% MAP estimate using all data
% initial hyperparameters
datastruct.muf = mean(r);
ovrscl_1 = max(r); % overall scale
lngthscl_1 = (max(max(support_x))-min(min(support_x)))*0.5; % variance
prs0 = [datastruct.muf; ovrscl_1; lngthscl_1];

[prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0] = updateFmapHyperparam_main(prs0, datastruct);

%% make prediction

[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
   
%%
load fmean_logexp1_cubic.mat
% load greg2D_alldata.mat
gfmean = fmean_logexp1_cubic(predictiveMean, sqrt(predictiveVar));
gfmean_2Dalldata = gfmean;
save gfmean_2Dalldata gfmean_2Dalldata

%% plot
% load gfmean_2Dalldata;
load  greg2D_alldata.mat
lengAxis = 100;
xaxis = linspace(min(support_x(:,1)), max(support_x(:,1)), lengAxis);
yaxis = linspace(min(support_x(:,2)), max(support_x(:,2)), lengAxis);

[xx, yy] = meshgrid(xaxis, yaxis);
% F = TriScatteredInterp(support_x(:,1), support_x(:,2), gfmean);
% qz = F(xx, yy);
% figure(2); contour(xx, yy, qz);
% 
%%
datastruct.support = [xx(:) yy(:)];

[predictiveMean_onGrid, predictiveVar_onGrid, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
qz = fmean_logexp1_cubic(predictiveMean_onGrid, sqrt(predictiveVar_onGrid));

figure(4); subplot(221); contour(xx, yy, reshape(qz, length(xx), []), 5); title('log-exp-g');
