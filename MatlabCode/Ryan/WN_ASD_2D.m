%% Implement 2D ASD on STA data for multiple frames and plot
%
% This function loads data for each offset value. For each offset value, it
% computes the isotropic and isotropic ASD estimates of the stimuli weights
% These are computed using 80% of the data. The OLS solution is also
% computed using the same 80% of the data. The weight outputs from the
% isotropic, anisotropic, and OLS are used to calculate the MSE on the held
% out 20% of the data. After calculating the MSE for each offset value, the
% results are plotted and saved if specified
%
%%
close all;
clear all;
clc;

%% Load data and define initialization parameters
% Load data
WN = nex2stro;
nstixperside = WN.sum.exptParams.nstixperside;
nframesidx = find(strcmp(WN.sum.trialFields(1,:), 'num_frames'));
nframes = sum(WN.trial(:,nframesidx));
basisvec = eye(3*nstixperside^2);
maxT = 1; % Number of frames to return. (More than 1 may exceed memory)

% Define path to save
p2s = ['D:\Users\Ryan\Desktop\' num2str(WN.sum.fileName(28:37)) ...
    '_LARS_Figs\'];

% Define offset values to analyze. 
offset = 0:8;

% ASD Parameters
nks = [10 10]; % Stimuli dimensions
nk = prod(nks); 
minlens = [1.5 1.5];  % minimum length scale along each dimension

% Center and standardize input and response if intended to compare to LARS
center_standard = true;

% Define rng seed. Set to [] to randomize, 'default' or and int to pick.
rng_seed = 'default';

% Save data and plot
save_data = true;
save_plots = true; 

%% Run ASD and MSE comparison
for ii = (offset+1)
    % .mex function init args
    initargs = {basisvec,ii-1, nframes,[nstixperside^2, 3, maxT]};
    
    % Call .mex function wrapper to get stimuli and responses
    projWN = getWhtnsStats(WN,maxT,'STPROJmod',initargs);
    
    % Separate data into training and test data
    if ~isempty(rng_seed); rng(rng_seed); end
    cvidx = crossvalind('KFold',size(projWN{1},1), 5);
    cvtrainidx = cvidx~=5;
    cvtestidx = cvidx == 5;
    % Center data before ASD
    X = projWN{1}(cvtrainidx, :)-mean(projWN{1}(cvtrainidx, :),1);
    y = projWN{2}(cvtrainidx,:)-mean(projWN{2}(cvtrainidx,:));
    Xtest = projWN{1}(cvtestidx,:)-mean(projWN{1}(cvtestidx,:),1);
    ytest = projWN{2}(cvtestidx,:)-mean(projWN{2}(cvtestidx,:));
    
    % Center and standardize (to unit norm) data if MSE is to be compared 
    % to LARS. 
    if center_standard
        % Scale and standardize x to unit length training data.
        X = X./sqrt(sum(X.^2,1));
        
        % Scale and standardize x to unit length training data.
        Xtest = Xtest./sqrt(sum(Xtest.^2,1));
    end
    
    % Compute ASD for each color individually and combine into one vector
    for jj = 0:2
        idx_start = (jj*100)+1; % Will pick index 1, 101, and 201 
        idx_end = (jj+1)*100; % Will pick index 100, 200, 300
        
%         % Call 2D ASD Isotropic
        [kasd(idx_start:idx_end,ii), asdstats] = ...
            fastASD(X(:,idx_start:idx_end), y, nks, minlens);
%         
%         kaldS(idx_start:idx_end,ii) = khatALD.khatS;
%         Call 2D ASD Anisotropic
        [kasd_aniso(idx_start:idx_end,ii),asdstats] = ...
            fastASD_aniso(X(:,idx_start:idx_end),y,nks,minlens);
    end
     
    % Calculate MSE for each frame isotropic and anisotropic
    CVmse_iso(ii) = sum((ytest-Xtest*kasd(:,ii)).^2,1)/size(Xtest, 1);
    CVmse_aniso(ii) = sum((ytest-Xtest*kasd_aniso(:,ii)).^2,1)/size(Xtest, 1);

    % Compute OLS MSE
    bols_train(:,ii) = X\y; % OLS from training data
    OLSmse(ii) = sum((ytest-Xtest*bols_train(:,ii)).^2,1)/size(Xtest, 1);
end

%% Plot Isotropic
WN_plotSTAest(kasd, nstixperside, 'ASD Iso Indv', save_plots, p2s);

%% Plot and scale Anisotropic results
WN_plotSTAest(kasd_aniso, nstixperside, 'ASD Aniso Indv', save_plots, p2s);

%% Plot MSE results
figure; hold on; 
plot(0:8,OLSmse,'LineWidth',2);
plot(0:8,CVmse_iso, 'LineWidth', 2); 
plot(0:8,CVmse_aniso,'LineWidth',2); 
hold off; 
legend('OLS STA', 'ASD Iso', 'ASD Aniso'); 
title(num2str(WN.sum.fileName(28:37)))
ylabel('MSE'); xlabel('Frame Offset');
if save_plots; saveas(gcf, [p2s 'LARS_STA_MSEcomp2.png']); end

% Save Data
if save_data
    save([p2s 'ASD_MSE.mat'], 'CVmse_iso', 'CVmse_aniso', 'OLSmse',...
        'offset', 'bols_train', 'kasd', 'kasd_aniso', 'center_standard');
end