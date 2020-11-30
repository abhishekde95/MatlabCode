% To find out which is the best space? Is it the space where the stimulus
% distribution is white. This is something we want to answer using this
% program
clear all;
close all;
%% Copied as it is from Greg's WNAnalysis
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

%% Gun Space - White Space
prompt1 = 'Do u want the stimulus distribution to be white in RG space or LM space? (1,2)\n';
str1 = input(prompt1,'s');
if (strcmp('1',str1))
   flag = 0;
elseif (strcmp('2',str1))
    flag =1 ;
end

mu = [0 0];
sigma = [1 0; 0 1];
num_points = 10000;
if (flag == 0)
    % creating points in the RG space using a bivariate normal distribution
    % points_RG is 10000 x 2 matrix
    % points_LM is 2 x 10000 matrix
    points_RG = mvnrnd(mu,sigma,num_points);
    points_LM = M(1:2,1:2)* points_RG'; % converting the points in RG dimension to LM dimension; (from gun space to cone space)
    contrast = [1 0] * points_LM; % projection onto the L+M line
    params = [50, 0.5, 2 1]; % A, sigmoid, exp, baseline
    FITSTR = 'symmetric';
    responses = ComputeNakaRushtonJPW(params,contrast,FITSTR); % This function is present in Patrick's directory
    zero_resp_ind = find(responses <= 0); % Find the indices of the stimuli which give no response
    STA_LM = sum(repmat(responses,[2 1]).* points_LM,2)/sum(responses);  % Computing the STA in LM space
    STA_RG = sum(repmat(responses',[1 2]).* points_RG,1)/sum(responses); % Computing the STA in RG space similar to the step above
    STA_RG_2_LM = M(1:2,1:2)*STA_RG'; % Converting the STA in LM domain to RG domain. This is done to check if STA_LM is same as STA_RG_2_LM
    
    % First subplot deals with points in RG space
    set(figure(1), 'Position', [100, 250, 1200, 450]);
    subplot(121);scatter(points_RG(:,1),points_RG(:,2),'g'); hold on;
    scatter(points_RG(zero_resp_ind,1),points_RG(zero_resp_ind,2),'r');
    scatter([0 STA_RG(1)],[0 STA_RG(2)]);
    line([0 STA_RG(1)],[0 STA_RG(2)],'LineWidth',2,'Color',[0 0 1.0]);
    line([-1*STA_RG(2) STA_RG(2)],[STA_RG(1) -1*STA_RG(1)],'LineWidth',2,'Color',[0 0 1.0]); % product of the slope of orthogonal lines are perpendicular
    xlabel('R (Red)'); ylabel('G (Green)'); title('Gun Space'); axis([-4 4 -4 4]); hold off;
    
    % Second subplot deals with points in LM subspace
    subplot(122); scatter(points_LM(1,:),points_LM(2,:),'g'); hold on;
    scatter(points_LM(1,zero_resp_ind),points_LM(2,zero_resp_ind),'r');
    scatter([0 STA_LM(1)],[0 STA_LM(2)]);
    line([0 STA_LM(1)],[0 STA_LM(2)],'LineWidth',2,'Color',[0 0 1.0]);
    line([0 STA_RG_2_LM(1)],[0 STA_RG_2_LM(2)],'LineWidth',2,'Color',[0 0 1.0]);
    line([-1*STA_LM(2) STA_LM(2)],[STA_LM(1) -1*STA_LM(1)],'LineWidth',2,'Color',[0 0 1.0]);
    xlabel('L (Long)'); ylabel('M (Medium)'); title('Cone Space'); axis([-0.75 0.75 -0.75 0.75]); hold off
else
    % creating points in the LM space using a bivariate normal distribution
    % points_LM is 10000 x 2 matrix
    % points_RG is 2 x 10000 matrix
    points_LM = mvnrnd(mu,sigma,num_points);
    points_RG = inv(M(1:2,1:2))* points_LM'; % converting the points in RG dimension to LM dimension; (from gun space to cone space)
    contrast = [1 -1] * points_RG; % projection onto the L+M line
    params = [50, 0.5, 2 1]; % A, sigmoid, exp, baseline
    FITSTR = 'symmetric';
    responses = ComputeNakaRushtonJPW(params,contrast,FITSTR); % This function is present in Patrick's directory
    zero_resp_ind = find(responses <= 0); % Find the indices of the stimuli which give no response
    STA_RG = sum(repmat(responses,[2 1]).* points_RG,2)/sum(responses);  % Computing the STA in LM space
    STA_LM = sum(repmat(responses',[1 2]).* points_LM,1)/sum(responses); % Computing the STA in RG space similar to the step above
    STA_LM_2_RG = inv(M(1:2,1:2))*STA_LM'; % Converting the STA in LM domain to RG domain. This is done to check if STA_LM is same as STA_RG_2_LM
    points_RG = points_RG';
    points_LM = points_LM';
    
    % First subplot deals with points in RG space
    set(figure(1), 'Position', [100, 250, 1200, 450]);
    subplot(121);scatter(points_RG(:,1),points_RG(:,2),'g'); hold on;
    scatter(points_RG(zero_resp_ind,1),points_RG(zero_resp_ind,2),'r');
    scatter([0 STA_RG(1)],[0 STA_RG(2)]);
    line([0 STA_RG(1)],[0 STA_RG(2)],'LineWidth',2,'Color',[0 0 1.0]);
     line([0 STA_LM_2_RG(1)],[0 STA_LM_2_RG(2)],'LineWidth',2,'Color',[0 0 1.0]);
    line([-1*STA_RG(2) STA_RG(2)],[STA_RG(1) -1*STA_RG(1)],'LineWidth',2,'Color',[0 0 1.0]); % product of the slope of orthogonal lines are perpendicular
    xlabel('R (Red)'); ylabel('G (Green)'); title('Gun Space'); axis([-100 100 -100 100]); hold off;
    
    % Second subplot deals with points in LM subspace
    subplot(122); scatter(points_LM(1,:),points_LM(2,:),'g'); hold on;
    scatter(points_LM(1,zero_resp_ind),points_LM(2,zero_resp_ind),'r');
    scatter([0 STA_LM(1)],[0 STA_LM(2)]);
    line([0 STA_LM(1)],[0 STA_LM(2)],'LineWidth',2,'Color',[0 0 1.0]);
    line([-1*STA_LM(2) STA_LM(2)],[STA_LM(1) -1*STA_LM(1)],'LineWidth',2,'Color',[0 0 1.0]);
    xlabel('L (Long)'); ylabel('M (Medium)'); title('Cone Space'); axis([-4 4 -4 4]); hold off
end

% Although I have calculated these variables but not explicitly used anywhere in the program
cov_RG = cov(points_RG(:,1),points_RG(:,2)); % calculates the variance and covariance between the R and the G samples
cov_LM = cov(points_LM(1,:),points_LM(2,:));
% The 2 points STA_LM and STA_RG_2_LM seem to overlap from the 2 plots