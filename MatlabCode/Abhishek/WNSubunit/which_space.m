% To find out which is the best space? Is it the space where the stimulus
% distribution is white.

clear all;
close all;
mu = [0 0];
sigma = [1 0; 0 1];
% creating points in the RG space;
points_RG = [];
num_points = 10000;
for i = 1:num_points
    points_RG = [points_RG; mvnrnd(mu,sigma)];
end
% display the points
figure(1);scatter(points_RG(:,1),points_RG(:,2)); xlabel('R'); ylabel('G');
cov_RG = cov(points_RG(:,1),points_RG(:,2)); % calculates the variance and covariance between the R and the G samples

origin = [0 0]; STA_RG = [0 0];
while(sum(abs(STA_RG)>0.5))
    STA_ind = randi(num_points);
    STA_RG = points_RG(STA_ind,:);
end
slope = STA_RG(2)/STA_RG(1);
slope_orth = -1/slope;

slope_RG_points = points_RG(:,2)./points_RG(:,1);
slopes_ind = find(abs(slope_RG_points - slope_orth)<0.01);
[~,max_ind] = max(points_RG(slopes_ind,1));
[~,min_ind] = min(points_RG(slopes_ind,1));
figure(2), scatter([origin(1) STA_RG(1)],[origin(2) STA_RG(2)]); hold on;
line([origin(1) STA_RG(1)],[origin(2) STA_RG(2)]);
scatter(points_RG(slopes_ind,1),points_RG(slopes_ind,2),'r'); xlabel('R'); ylabel('G'); 
line([points_RG(slopes_ind(min_ind),1) points_RG(slopes_ind(max_ind),1)],[points_RG(slopes_ind(min_ind),2) points_RG(slopes_ind(max_ind),2)]);
axis([-3 3 -3 3]);
hold off;

RG_transform = [origin(1) STA_RG(1) points_RG(slopes_ind,1)';origin(2) STA_RG(2) points_RG(slopes_ind,2)' ];

%%
% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
stro = nex2stro(findfile('N021915003.nex'));
fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');

%% 
% Analysing in LM space where the distribution of the stimuli is not white 
points_LM_all = M(1:2,1:2) * points_RG';
points_LM = M(1:2,1:2)* RG_transform;
[~,maxLM_ind] = max(points_LM(1,3:end));
[~,minLM_ind] = min(points_LM(1,3:end));
figure(3),scatter(points_LM_all(1,:), points_LM_all(2,:),'g'); hold on;
scatter(points_LM(1,1:2), points_LM(2,1:2)); line(points_LM(1,1:2), points_LM(2,1:2));
scatter(points_LM(1,3:end),points_LM(2,3:end),'r'); xlabel('L'); ylabel('M'); hold off; 
