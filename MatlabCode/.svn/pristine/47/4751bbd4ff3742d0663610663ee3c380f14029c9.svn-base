%% Run LARS and Plot

%% Diabetes Data from Efron et al. 2004
data = readmatrix('D:\Users\Ryan\Desktop\LARS_diabetes_dataset.csv');
[beta_efron, info_efron] = WN_LARS(data(:,1:10), data(:,11), true, true);

% Plot weights (match fig 3-left) from paper
figure; plot(sum(abs(beta_efron),1),beta_efron', 'linewidth', 3), 
lgd = legend('1', '2', '3', '4','5','6','7','8','9','10', 'Location', ...
    'westoutside');
lgd.FontSize = 14;
% xlim([0, 175])
title('LARS: Data from Efron et al.', 'FontSize', 18);
xlabel('t = sum(|\beta_j|)', 'FontSize', 16); ylabel('\beta_j', ...
    'FontSize', 16);

% Plot Cp (match fig 7-left) from paper
figure; plot(info_efron.cp,'k', 'linewidth', 3);
xlim([4.8, 10.2]); 
xticks([5, 6, 7, 8, 9, 10]);
title('LARS Cp: Data from Efron et al.', 'FontSize', 18);
xlabel('Step Number', 'FontSize', 16); ylabel('Estimated Cp', 'FontSize', 16);

% Compare to OLS solution
beta_efronOLS = beta_efron(:,end)./max(beta_efron(:,end));
%preprocess
rawstim = data(:,1:10) - mean(data(:,1:10));
for kk = 1:size(data(:,1:10),2)
    efronstim(:,kk) = rawstim(:,kk)/norm(rawstim(:,kk));
end

efronresp = data(:,11) - mean(data(:,11));

% OLS
efronOLS = efronstim\efronresp;
efronOLS = efronOLS./max(efronOLS);

%% Simulated WN Data
nstim = 300;
n = 10;
rng(150)
stimuli = normrnd(0,1,nstim,n);

% Simulated responses
rng(5)
resp_act = normrnd(0,3,1,nstim);

% Call LARS
tic
[beta_simWN, info_simWN] = WN_LARS(stimuli, resp_act, true);
toc

% Plot weights
figure; plot(sum(abs(beta_simWN),1),beta_simWN', 'linewidth', 3), 
lgd = legend('1', '2', '3', '4','5','6','7','8','9','10', 'Location', 'westoutside');
lgd.FontSize = 14;
xlim([0, 2.1])
title('LARS: Simulated Data', 'FontSize', 18);
xlabel('t = sum(|\beta_j|)', 'FontSize', 16); ylabel('\beta_j', 'FontSize', 16);

% Plot Cp 
figure; plot(info_simWN.cp,'k', 'linewidth', 3); 
xticks([1:length(info_simWN.cp)]);
title('LARS Cp: Simulated Data', 'FontSize', 18);
xlabel('Step Number', 'FontSize', 16); ylabel('Estimated Cp', 'FontSize', 16);

%% Processing time test both:
nstim = 500;
dims2test = 5:10:490;
fminsch_time = zeros(length(dims2test),1);
lars_time = zeros(length(dims2test), 1);

for nn = 1:length(dims2test)
    % Define X and y
    rng(12)
    stimuli = normrnd(0,1,nstim,dims2test(nn));
    rng(7)
    resp_act = normrnd(0,3,nstim, 1);
    
    % Call fminsearch
    tic 
    [~] = weights_test(stimuli, resp_act, true);
    fminsch_time(nn) = toc;
    
    % Call LARS
    tic
    [~, ~] = WN_LARS(stimuli, resp_act, true);
    lars_time(nn) = toc;
    disp(num2str(100*nn/length(dims2test)));
    % Display percent complete
    if mod(nn,5) == 0
        disp(num2str(100*nn/length(dims2test)));
    end
end

% Plot
figure; 
semilogy(dims2test, fminsch_time, 'LineWidth', 4);hold on;
semilogy(dims2test, lars_time, 'LineWidth', 4);
xlabel('Dimensions', 'FontSize', 16); ylabel('Elapsed Time [s]', 'FontSize', 16);
legend('fminsearch', 'LARS', 'Location', 'northwest','FontSize', 14);
title('Elapsed Time comparison, nstim = 500', 'FontSize', 18);
hold off; 

%% Processing time test LARS:
nstim = 500;
dims2test = 5:10:490;
fminsch_time = zeros(length(dims2test),1);
lars_time2 = zeros(length(dims2test), 1);

for nn = 1:length(dims2test)
    % Define X and y
    rng(12)
    stimuli = normrnd(0,1,nstim,dims2test(nn));
    rng(7)
    resp_act = normrnd(0,3,nstim, 1);
    
    % Call LARS
    tic
    [~, ~] = WN_LARS(stimuli, resp_act, true);
    lars_time2(nn) = toc;
    disp(num2str(100*nn/length(dims2test)));
    % Display percent complete
    if mod(nn,5) == 0
        disp(num2str(100*nn/length(dims2test)));
    end
end

% Plot
figure; 
semilogy(dims2test, lars_time, 'LineWidth', 4);
xlabel('Dimensions', 'FontSize', 16); ylabel('Elapsed Time [s]', 'FontSize', 16);
legend('LARS', 'Location', 'northwest','FontSize', 14);
title('Elapsed Time, nstim = 500', 'FontSize', 18);
hold off; 
