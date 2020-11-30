function [next_stims,pred_threshs] = next_directions_to_test(x1, x2, m, ps2, nstim, hyp, x, y)
%  [next_stims,pred_threshs] = next_directions_to_test(x1, x2, m, ps2, nstim)
%
%  Selects up to nstim new directions to probe in an LMTF experiment.
%  Stimulus direction selection is based uncertainty sampling using a
%  Gaussian process regression model.
%
%  INPUT (all of these inputs except nstim come from evaluate_posterior.m)
%      x1: a meshgrid of log10(temporal frequencies) from log10(minTF) to
%      log10(maxTF)
%      x2: a meshgrid of color direcitons from 0 to pi 
%      m: the mean of the GP fit at each (x1, x2) pair
%      ps2: the variance of the fit at each (x1, x2) pair
%      nstim: the maximum number of stimuli to select
%
%  OUPUT
%      next_stims: nx2 matrix of (TF, theta) pairs to be tested.
%      pred_threshs nx1 vector of predicted thresholds for each of the
%      newly selected stimulus directions.
%
%   Algorithm:
%      Select the point on the 50x50 (TFxtheta) lattice at which ps2 is
%      greatest (the point at which our estimate of the mean function is
%      worst). Then, pin that point at the mean function, re-evaluated the
%      posterior, and recurse.
%
%   GDLH 5/11/16


%   Old algorithm:
%      We sort the ps2s in descending order to find the point at which 
%      where ps2 is greatest (our estimate of the mean function is worst).
%      That point is selected. Then we go down the list, never taking a
%      point that neighbors a point that we've already considered. This
%      biases us towards taking points that are local maxima in ps2
%      (imagine dropping a horzontal plane onto the ps2 surface; every time
%      a new spot appears in the intersection, we select that stimulus
%      condition for testing.) This means we tend to sample from different
%      "hills" in in the ps2 surface, and we will not sample along ridges
%      of high ps2. A reasonable alternative would be to sample from the
%      model and recurse.
%

%NBHD_RADIUS = 1; % NBHD = "neighborhood"  % Needed for old algorithm
[~,scalars] = gamut_extent([cos(x2(:,1)) sin(x2(:,1))]);
gamutmask = repmat(log10(scalars'), 1, size(x2, 2));
predOOGpoints = m>gamutmask(:); % points that are predicted to be out of gamut
ps2(predOOGpoints) = 0; % if the GP fit is outside of the gamut, don't even try to test there
next_stims = zeros(nstim,2);
pred_threshs = zeros(nstim,1);
mask_size = size(gamutmask); % 50 x 50 

% OLD METHOD FOR PICKING STIMULI
% ------------------------------
%
% mask = zeros(mask_size);
% pred_threshs = zeros(nstim,1);
% var_idx = 1;
% stim_idx = 1;
% 
% [~,I] = sort(ps2, 'descend');
% while stim_idx <= nstim && var_idx < length(ps2)
%     [ii,jj] = ind2sub(mask_size, I(var_idx));
%     nbhd_rows = mod((-NBHD_RADIUS:NBHD_RADIUS)+ii-1, mask_size(1))+1; % the rows are periodic
%     nbhd_cols = (-NBHD_RADIUS:NBHD_RADIUS)+jj;
%     nbhd_cols(nbhd_cols < 1) = []; % Don't want to index ouside the range
%     nbhd_cols(nbhd_cols > mask_size(2)) = []; % Don't want to index ouside the range
%     if ~any(any(mask(nbhd_rows,nbhd_cols))) % if we're not within NBHD_RADIUS of a point already selected this round
%         next_stims(stim_idx,:) = [x1(ii,jj) x2(ii,jj)];
%         pred_threshs(stim_idx) = m(ii+(jj-1)*mask_size(1));
%         stim_idx = stim_idx+1;
%     end
%     mask(ii,jj) = 1;
%     var_idx = var_idx+1;
% end
% next_stims(stim_idx:end,:) = [];
% pred_threshs(stim_idx:end) = [];


% NEW METHOD FOR PICKING STIMULI
% ------------------------------
minTF = x1(1,1);
maxTF = x1(1,end);
for stim_idx = 1:nstim 
    % For testingminTF
    %figure; 
    %subplot(2,1,1); surface(10.^x1,x2,reshape(m,mask_size));
    %subplot(2,1,2); surface(10.^x1,x2,reshape(ps2,mask_size));
    % End of testing code
    whichpoint = find(ps2 == max(ps2));
    [ii,jj] = ind2sub(mask_size,whichpoint(unidrnd(length(whichpoint))));
    next_stims(stim_idx,:) = [x1(ii,jj) x2(ii,jj)];
    pred_threshs(stim_idx) = m(ii+(jj-1)*mask_size(1));
    x = [x; next_stims(stim_idx,:)];
    y = [y; pred_threshs(stim_idx)];
    [~, ~, m, ps2] = evaluate_posterior(hyp, x, y, minTF, maxTF);
    ps2(predOOGpoints) = 0;
end