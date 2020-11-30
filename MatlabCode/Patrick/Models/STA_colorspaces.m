% To find out which is the best space? Is it the space where the stimulus
% distribution is white. This is something we want to answer using this
% program
clear all;
close all;
%% Copied as it is from Greg's WNAnalysis
% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
%stro = nex2stro(findfile('N021915003.nex'));

% JPW edit to find file
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
end
datafile = 'M061515001';
stro = nex2stro([library char(datafile) '.nex']);

fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
%mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = stro.sum.exptParams.monspd; % JPW edit
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = M(1:2,1:2); % JPW edit

%% Gun Space - White Space
% creating points in the RG space using a bivariate normal distribution
% points_RG is 10000 x 2 matrix
% points_LM is 2 x 10000 matrix
mu = [0 0];
sigma = [1 0; 0 1];
%points_RG = [];
num_points = 10000;
%for i = 1:num_points
%    points_RG = [points_RG; mvnrnd(mu,sigma)];
%end
points_RG = mvnrnd(mu,sigma,num_points); % vectorizing to avoid loop (JPW)
%points_LM = M(1:2,1:2)* points_RG'; % converting the points in RG dimension to LM dimension; (from gun space to cone space)
points_LM = points_RG * M'; % JPW edit


% Haave calculated these variables but have not explicitly used anywhere in the program
cov_RG = cov(points_RG(:,1),points_RG(:,2)); % calculates the variance and covariance between the R and the G samples
cov_LM = cov(points_LM(:,1),points_LM(:,2));

artificial_tuning = [1 0 i mean, i]; % explicitly naming this variable (JPW)
%contrast = [1 1] * points_LM; % projection onto the L+M line
contrast = [artificial_tuning * points_LM']'; % JPW edit
params = [50, 0.5, 2 1]; % A, sigmoid, exp, baseline
FITSTR = 'symmetric';
responses = ComputeNakaRushtonJPW(params,contrast,FITSTR); % This function is present in Patrick's directory
%responses = contrast;
zero_resp_L = (responses <= 0); % Find the indices of the stimuli which give no response 

%STA_LM = mean(repmat(responses,[2 1]).* points_LM,2)/mean(responses);  % Computing the STA in LM space
STA_LM = sum(points_LM .* repmat(responses,1,2),1)./sum(responses)*2; % JPW edit
%STA_RG = mean(repmat(responses',[1 2]).* points_RG,1)/mean(responses); % Computing the STA in RG space similar to the step above
STA_RG = sum(points_RG .* repmat(responses,1,2),1)./sum(responses)*2; % JPW edit

%STA_RG_2_LM = M(1:2,1:2)*STA_RG; % Converting the STA in LM domain to RG domain. This is done to check if STA_LM is same as STA_RG_2_LM
STA_RG_2_LM = inv(M')*STA_RG'; % Greg's transform to preserve the mechanism

% Visualize the results - Display the points
figure(1); clf;
set(gcf,'Pos',[100 250 1200 600]);

% Stimuli and STA in gun space
axes('parent',gcf,'units','normalized','pos',[.2 .4 .25 .5],'box','on');
hold on; axis equal; grid on;
scatter(points_RG(~zero_resp_L,1),points_RG(~zero_resp_L,2),'g'); % JPW edit
scatter(points_RG(zero_resp_L,1),points_RG(zero_resp_L,2),'r');
plot([0 STA_RG(1)],[0 STA_RG(2)],'b'); 
plot([-1*STA_RG(2) STA_RG(2)],[STA_RG(1) -1*STA_RG(1)],'b'); % product of the slope of orthogonal lines are perpendicular
lim = max(abs(points_RG(:)));
xlim([-lim lim]); ylim([-lim lim])
xlabel('R (Red)'); ylabel('G (Green)'); title('Gun Space');

% Red gun histogram
axes('parent',gcf,'units','normalized','pos',[.2 .1 .25 .2],'box','on');
[counts,centers] = hist(points_RG(:,1),50); 
bar(centers,counts); hold on;
plot([mean(points_RG(:,1)) mean(points_RG(:,1))],[0 max(counts)],'r-')
xlim([-lim lim])

% Green gun histogram
axes('parent',gcf,'units','normalized','pos',[.05 .4 .1 .5],'box','on');
[counts,centers] = hist(points_RG(:,2),50);
barh(centers,counts); hold on;
plot([0 max(counts)],[mean(points_RG(:,2)) mean(points_RG(:,2))],'r-')
ylim([-lim lim])

% Stimuli and STA in cone space
axes('parent',gcf,'units','normalize','pos',[.5 .4 .25 .5],'box','on');
hold on; axis equal; grid on;
scatter(points_LM(~zero_resp_L,1),points_LM(~zero_resp_L,2),'g');
scatter(points_LM(zero_resp_L,1),points_LM(zero_resp_L,2),'r');
plot([0 STA_LM(1)],[0 STA_LM(2)],'b');
plot([-1*STA_LM(2) STA_LM(2)],[STA_LM(1) -1*STA_LM(1)],'b');
plot([0 STA_RG_2_LM(1)],[0 STA_RG_2_LM(2)],'color',[1 1 0]);
xlabel('L (Long)'); ylabel('M (Medium)'); title('Cone Space'); 
xlim([-.6 .6]); ylim([-.6 .6])

% L-cone histogram
axes('parent',gcf,'units','normalized','pos',[.5 .1 .25 .2]);
[counts,centers] = hist(points_LM(:,1),50); 
bar(centers,counts); hold on;
plot([mean(points_LM(:,1)) mean(points_LM(:,1))],[0 max(counts)],'r-')
xlim([-.6 .6])

% M-cone histogram
axes('parent',gcf,'units','normalized','pos',[.8 .4 .1 .5]);
[counts,centers] = hist(points_LM(:,2),50);
barh(centers,counts); hold on;
plot([0 max(counts)],[mean(points_LM(:,2)) mean(points_LM(:,2))],'r-')
ylim([-.6 .6])

% The 2 points STA_LM and STA_RG_2_LM seem to overlap from the 2 plots