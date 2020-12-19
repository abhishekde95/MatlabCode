%% WN RUN ALD

%% Load data

% Load data
WN = nex2stro;
nstixperside = WN.sum.exptParams.nstixperside;
nframesidx = find(strcmp(WN.sum.trialFields(1,:), 'num_frames'));
nframes = sum(WN.trial(:,nframesidx));
basisvec = eye(3*nstixperside^2);
maxT = 1; % Number of frames to return from STPROJmod. (More than 1 may exceed memory)

% Define path to save
exp_name = num2str(WN.sum.fileName(28:37));
p2s = ['D:\Users\Ryan\Desktop\' exp_name '_LARS_Figs\'];

% Define offset values to analyze. 
offset = 0:8;

% Center and standardize input and response if intended to compare to LARS
cross_val = true;
center_standard = false;

% Define rng seed. Set to [] to randomize, 'default' or and int to pick.
rng_seed = 'default';

% Save data and plot
save_data = true;
save_plots = true; 

zerosidx = zeros(length(offset), 3);

%% Run ALD and MSE comparison
for ii = (offset+1)
    % .mex function init args
    initargs = {basisvec,ii-1, nframes,[nstixperside^2, 3, maxT]};
    
    % Call .mex function wrapper to get stimuli and responses
    projWN = getWhtnsStats(WN,maxT,'STPROJmod',initargs);
    
    if cross_val
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
    else
        X = projWN{1};
        y = projWN{2};
    end
    
    % Compute ASD for each color individually and combine into one vector
    for jj = 0:2
        idx_start = (jj*100)+1; % Will pick index 1, 101, and 201 
        idx_end = (jj+1)*100; % Will pick index 100, 200, 300
         
        [khatALD, kRidge]= runALD(X(:,idx_start:idx_end),...
            y, [nstixperside;nstixperside], 1);
        
        ALDs(idx_start:idx_end,ii) = khatALD.khatS;
%         ALDf(idx_start:idx_end,ii) = khatALD.khatF;
        ALDsf(idx_start:idx_end,ii) = khatALD.khatSF;
        
        % Track which stimuli/responses cause kRidge = 0
        if all(kRidge == 0)
            disp(['kRidge = 0: Offset ' num2str(ii-1) ', Color ' num2str(jj)])
            zerosidx(ii+1, jj+1) = 1;
        end        
    end
    
    if cross_val
        % Calculate MSE for each frame isotropic and anisotropic
        CVmse_ALDs(ii) = sum((ytest-Xtest*ALDs(:,ii)).^2,1)/size(Xtest, 1);
        CVmse_ALDsf(ii) = sum((ytest-Xtest*ALDsf(:,ii)).^2,1)/size(Xtest, 1);
        
        % Compute OLS MSE
        bols_train(:,ii) = X\y; % OLS from training data
        OLSmse(ii) = sum((ytest-Xtest*bols_train(:,ii)).^2,1)/size(Xtest, 1);
    end
end


%% Plot
WN_plotSTAest(ALDs, nstixperside, 'ALDs', save_plots, p2s)
WN_plotSTAest(ALDsf, nstixperside, 'ALDsf', save_plots, p2s)


% Plot MSE
if cross_val
    figure; hold on;
    plot(offset, CVmse_ALDs, 'LineWidth', 3);
    plot(offset, CVmse_ALDsf, 'LineWidth', 3);
    plot(offset, OLSmse, 'LineWidth', 3);
    legend('ALDs', 'ALDsf', 'OLS');
    title('CV MSE Comp');
    xlabel('Offset');
    ylabel('MSE');
    hold off;
    
    if save_plots; saveas(gcf, [p2s 'CV_MSE_ALD.png']); end
end
% Save Data
if save_data && cross_val
    save([p2s 'ALD_MSE.mat'], 'ALDs', 'ALDsf', 'offset', 'exp_name',...
        'CVmse_ALDs', 'CVmse_ALDsf', 'OLSmse', 'zerosidx');
    
elseif save_data
    save([p2s 'ALD.mat'], 'ALDs', 'ALDsf', 'offset', 'exp_name', 'zerosidx')
end