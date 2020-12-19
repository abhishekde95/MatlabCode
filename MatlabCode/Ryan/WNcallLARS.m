%% White noise LARS
% This script calls LARS for each frame offset and plots the Cp value at
% each offset, the resulting weights, and the OLS solution. 
%%

close all
clear all
clc

%% Load WN data
WN=nex2stro;

% Path to save dataa
p2s = ['D:\Users\Ryan\Desktop\' num2str(WN.sum.fileName(28:37)) ...
    '_LARS_Figs\'];
save_figs = false;

% LARS Steps to plot
lars_steps = [51, 101, 151, 301];

% Run and plot LARS for a range of offset values.
nstixperside = WN.sum.exptParams.nstixperside;
nframesidx = find(strcmp(WN.sum.trialFields(1,:), 'num_frames'));
nframes = sum(WN.trial(:,nframesidx));
basisvec = eye(3*nstixperside^2);
maxT = 1;

% Frames 
frame_offset = 0:8;

%% Run LARS at each offset and plot Cp for each offset
figure('Renderer', 'painters', 'Position', [100 100 2000 200])
for ii = (frame_offset+1)
    
    % .mex function init args
    initargs = {basisvec,ii, nframes,[nstixperside^2, 3, maxT]};
    
    % Call .mex function wrapper to get stimuli and responses
    projWN = getWhtnsStats(WN,maxT,'STPROJmod',initargs);
   
    % Preprocess raw data to run regression control
    prestim = projWN{1}-mean(projWN{1});
    for kk = 1:size(projWN{1},2)
        prestim(:,kk) = prestim(:,kk)/norm(prestim(:,kk));
    end

    % Run regression control
    projWNregress(:,ii) = prestim\(projWN{2}-mean(projWN{2}));
    clear prestim
    
    % Run LARS on data
    tic
    [beta, fitInfo] = WN_LARS(projWN{1}, projWN{2}, true, true);
    toc
    
    % Clear 
    clear projWN prestim
    % Save data we're interested in
    betaAll(:,:,ii) = beta;
    cp(:,ii) = fitInfo.cp;
    
    subplot(1,9,ii)
    plot(fitInfo.cp, 'k', 'LineWidth', 2);
    title(['LARS Cp - Offset: ' num2str(ii)], 'FontSize', 12);
    xlabel('Step', 'FontSize', 10);
    box off;
end
if save_figs; saveas(gcf, [p2s 'LARS_Cp.png']); end
disp('Done')

%% Save the results at discrete LARS steps for all frames
for kk = 1:length(lars_steps)
    WN_plotSTAest(squeeze(betaAll(:,lars_steps(kk),:)), ...
        nstixperside, ['LARS ' num2str(lars_steps(kk))], save_figs, p2s)
end

%% Plot Regrssion Control
WN_plotSTAest(projWNregress,nstixperside,'STPROJmod OLS', save_figs, p2s)

%% Calcalate and plot each offset based on the minCp value
[betaParse, minCpidx] = LARS_minCp(betaAll, cp);

% Scale and plot
WN_plotSTAest(betaParse,nstixperside,'STPROJmod OLS', save_figs, p2s)

%% Preprocess raw data
% rawstim = projWN{1}-mean(projWN{1});
% for kk = 1:size(projWN{1},2)
%     prestim(:,kk) = rawstim(:,kk)/norm(rawstim(:,kk));
% end
% 
% preresp = projWN{2}-mean(projWN{2});