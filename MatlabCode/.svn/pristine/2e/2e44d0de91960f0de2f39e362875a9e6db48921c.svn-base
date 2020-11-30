% main function to test active learning algorithm using persistent variables

clear;
clc;

maxNumbTrials = 500; % number of datapoints we will observe by active learning
stimDim = 2;
numInitData = stimDim*10; % initial datapoints
totNumTrials = maxNumbTrials + numInitData -1;  % total number of trials

%% true lambda
npts = 25;
support_x = 1:npts;
support_y = 1:npts;
[xx, yy] = meshgrid(support_x, support_y);

support = [xx(:) yy(:)];

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

%% initial stimuli

numInitData = stimDim*10;

% initial points from LH design
stim_range_min = min(support);
stim_range_max = max(support);
x_lh = lhsdesign(numInitData, stimDim);
xInit = bsxfun(@plus, bsxfun(@times, (stim_range_max - stim_range_min), x_lh), stim_range_min);

rInit = poissrnd(g(intensity(xInit)));
rInit(rInit==0) = 0.5;

figure(300); clf;
subplot(221); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; title('true firing map');
subplot(222); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; hold on;
plot(xInit(:,1), xInit(:,2), 'o'); title('initial points from LH');

%%

count = numInitData;

% build a data structure (only once in its first trial after 20 initial
% observations) 
datastruct.x = xInit;
datastruct.r = rInit;
datastruct.support = support;
datastruct.ndim = stimDim;
datastruct.numInitData = numInitData; 
datastruct.totNumTrials = totNumTrials; 
datastruct.norm_mat_support = form_normMat(support, support);  % squared distance

g = @(t) log(exp(t)+1);
ginv = @(t) log(exp(t)-1);
datastruct.g = g;
datastruct.ginv = ginv;
datastruct.finit = datastruct.ginv(datastruct.r+0.1);

datastruct.thrsh_detH = 0.005;
datastruct.detH = 1;

datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);

load fvar_logexp1_lin.mat;
load fmean_logexp1_lin.mat;

datastruct.fvar_logexp1 = fvar_logexp1_lin;
datastruct.fmean_logexp1 = fmean_logexp1_lin;

% delete this for real experiments
datastruct.mse = zeros(totNumTrials-numInitData,1);
datastruct.lambdatrue = lambda_true;

% start active learning 

while count<=totNumTrials
    
    [ count totNumTrials ]
    
        [nextX, datastruct] = computeNextStim_ALalgorithm_datastruct(datastruct);
        rNew = poissrnd(g(intensity(nextX)));
        
        count = count +1;
        % update datastructure
        datastruct.x = [datastruct.x; nextX];
        datastruct.r = [datastruct.r; rNew+0.1]; % to avoid numerical problems, add some nonzero value to zero spike count.
       
        sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
        datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
        
        normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
        datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
        datastruct.finit = [datastruct.finit; datastruct.ginv(datastruct.r(end)+0.1)];        
        
end

