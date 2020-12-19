%% White noise LARS with cross validation for a single frame
% close all
% clear all
% clc

%% Load WN data
WN=nex2stro;

% Path to save dataa
p2s = ['D:\Users\Ryan\Desktop\' num2str(WN.sum.fileName(28:37)) ...
    '_LARS_Figs\'];
save_figs = true; 

% Run and plot LARS for a range of offset values.
nstixperside = WN.sum.exptParams.nstixperside;
nframesidx = find(strcmp(WN.sum.trialFields(1,:), 'num_frames'));
nframes = sum(WN.trial(:,nframesidx));
basisvec = eye(3*nstixperside^2);
maxT = 1;
kfold = 5;
frame_offset = 0:8;

% Initialize MSE results
minStep = zeros(length(frame_offset),1);

for kk = 1:length(frame_offset)
    % .mex function init args
    initargs = {basisvec,frame_offset(kk), nframes,[nstixperside^2, 3, maxT]};
    
    % Call .mex function wrapper to get stimuli and responses
    projWN = getWhtnsStats(WN,maxT,'STPROJmod',initargs);

    % K-fold cross validation
    [minStep(kk), mseStep(kk,:),~,cvidx] = WNLARS_kfoldcv(projWN{1}, ...
        projWN{2},kfold,[], true, true);

    % Run LARS on 80% of CV data to get model features
    [beta, ~] = WN_LARS(projWN{1}(cvidx~=kfold,:),...
        projWN{2}(cvidx~=kfold), true, true);
    
    % Center and normalize held out data
    X = projWN{1}(cvidx == kfold,:);
    y = projWN{2}(cvidx == kfold);
    Xc = X - mean(X,1);
    Xc = Xc./sqrt(sum(Xc.^2,1));
    yc = y-mean(y);
    
    % Extract weights for the 'best' step size given k-fold cv
    beta_cvopt(:,kk) = beta(:,minStep(kk));
    
    % Calculate MSE of CV features with 
    CVmse(kk) = sum((yc-Xc*beta(:,minStep(kk))).^2,1)/size(Xc, 1);
    
    % Calculate MSE of all features to compare
    STAmse(kk) = sum((yc-Xc*beta(:,end)).^2,1)/size(Xc, 1);
    
    disp(['Completed Lag ' num2str(kk)]);
end
disp('Done');


%% Plot MSE results
figure; hold on; 
plot(0:8,CVmse, 'LineWidth', 2); 
plot(0:8,STAmse,'LineWidth',2); 
hold off; 
legend('CV', 'STA'); 
title(num2str(WN.sum.fileName(28:37)))
ylabel('MSE'); xlabel('Frame Offset');
saveas(gcf, [p2s 'LARS_STA_MSEcomp.png']);

%% Save Data
save([p2s 'LARS_MSE.mat'], 'CVmse', 'STAmse', 'frame_offset', ...
    'beta_cvopt', 'minStep', 'mseStep', 'cvidx');

%% Plot beta at each 
WN_plotSTAest(beta_cvopt, nstixperside, 'Beta Optimal', save_figs, p2s)