% main function to test active learning algorithm using persistent variables

clear all;
clc;

maxNumbTrials = 500; % number of datapoints we will observe by active learning
stimDim = 2;
numInitData = stimDim*10; % initial datapoints
totNumTrials = maxNumbTrials + numInitData -1;  % total number of trials
x = zeros(totNumTrials, stimDim);
r = zeros(totNumTrials,1);

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

x(1:numInitData, :) = xInit;
r(1:numInitData) = poissrnd(intensity(xInit));
r(r==0) = 0.5;

figure(301); clf;
subplot(221); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; title('true firing map');
subplot(222); contour(xx, yy, reshape(lambda_true, npts, []), 5); axis image; axis xy; hold on;
plot(xInit(:,1), xInit(:,2), 'o'); title('initial points from LH');

%%

count = numInitData;
whichMethod = 1; % minimum variance
% whichMethod = 2; % max info

while count<=totNumTrials
    [ count totNumTrials ]
        nextX = computeNextStim_ALalgorithm(x(1:count, :), r(1:count), stimDim, support, numInitData, totNumTrials, lambda_true, whichMethod);
        rNew = poissrnd(intensity(nextX));
        
        count = count +1;
        x(count,:) = nextX;
        r(count) = rNew;
        r(r==0) = 0.5; % to avoid numerical problems, add some nonzero value to zero spike count.
       
end

