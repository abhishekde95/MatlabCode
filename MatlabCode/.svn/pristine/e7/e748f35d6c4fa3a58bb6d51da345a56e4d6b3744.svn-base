% For analyzing the cone spectral correlations
close all; clearvars;
load fundamentals.mat
load mon_spd.mat
cone_fundamentals = reshape(fundamentals,[numel(fundamentals)/3 3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = cone_fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
% M = inv(M');

% plotting the cone fundamentals and the scatter plots of pair of cone activations
figure(1),plot(linspace(380,780,numel(fundamentals)/3),cone_fundamentals,'Linewidth',2); xlabel('Wavelength'); ylabel('Power');


% converting RGB gaussian uncorrelated stimuli into LMS space 
mu = [0 0 0];
sigma = [1 0 0;0 1 0;0 0 1];
rnd_points = mvnrnd(mu,sigma,10000);
R = rnd_points(:,1);
G = rnd_points(:,2);
B = rnd_points(:,3);
RGB = cat(1,R',G',B');
LMS = M*RGB;
figure(2), subplot(221), scatter3(LMS(1,:),LMS(2,:),LMS(3,:)); xlabel('L'),ylabel('M'),zlabel('S');
subplot(222), scatter(LMS(1,:),LMS(2,:)); xlabel('L'),ylabel('M');
subplot(223), scatter(LMS(2,:),LMS(3,:)); xlabel('M'),ylabel('S');
subplot(224), scatter(LMS(1,:),LMS(3,:)); xlabel('L'),ylabel('S');

% Obtaining the correlations of the S/LM, L/M, Orange-cyan and Lime-magenta directions.
corr_LM = corr(LMS(1,:)',LMS(2,:)');
corr_OrangeCyan = corr(LMS(2,:)'-LMS(1,:)',LMS(3,:)'); % M+S-L
corr_LimeMagenta = corr(LMS(1,:)'-LMS(2,:)',LMS(3,:)'); % L+S-M
corr_S_LM = corr(-(LMS(1,:)'+LMS(2,:)'),LMS(3,:)'); % S-(L+M)
x = {'LM','M-L/S','L-M/S','S/LM'};
y = [corr_LM,corr_OrangeCyan,corr_LimeMagenta,corr_S_LM];
figure(3), bar(y); set(gca,'XTickLabel',x); ylabel('Stimulus Correlations'); title('WN RGB converted to LMS'); 

% finding out how much variance (of light) was delivered along different color directions 
lum_kernel = [1 1 0]; lum_vals = lum_kernel * LMS; 
rg_kernel = [1 -1 0]; rg_vals = rg_kernel * LMS;
Orange_Cyan_kernel = [1 -1 -1]; Orange_Cyan_vals = Orange_Cyan_kernel * LMS;
Lime_Magenta_kernel = [1 -1 1]; Lime_Magenta_vals = Lime_Magenta_kernel * LMS;
S_LM_kernel = [-1 -1 1]; S_LM_vals = S_LM_kernel * LMS;
figure(4), subplot(231); hist(lum_vals); xlabel('Contrast'); title('L+M'); text(0.4,2500, strcat('var=',num2str(var(lum_vals))));
subplot(232); hist(rg_vals); xlabel('Contrast'); title('L-M'); text(0.05,2500, strcat('var=',num2str(var(rg_vals))));
subplot(233); hist(Orange_Cyan_vals); xlabel('Contrast'); title('L-M-S (Orange-Cyan)'); text(0.2,2500, strcat('var=',num2str(var(Orange_Cyan_vals))));
subplot(234); hist(Lime_Magenta_vals); xlabel('Contrast'); title('L-M+S (Lime-Magenta)'); text(0.2,2500, strcat('var=',num2str(var(Lime_Magenta_vals))));
subplot(235); hist(S_LM_vals); xlabel('Contrast'); title('S-(L+M)'); text(0.3,2500, strcat('var=',num2str(var(S_LM_vals))));
subplot(236); scatter(LMS(1,:) - LMS(2,:), LMS(3,:)); xlabel('L-M'); ylabel('S'); title('Isoluminant plane');
