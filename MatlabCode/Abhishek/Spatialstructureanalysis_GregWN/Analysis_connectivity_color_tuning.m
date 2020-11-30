% Analysis of connectivity between LGN and V1, and color tuning of V1 cells
% Author - Abhishek De, 4/20

% Section 1: Figuring out the approximate number of LGN cells

% Section 2: Simulating resulting color tuning with random L- and M-weighted input

%% Section 1: Figuring out the approximate number of LGN cells 

close all; clearvars;
plot_counter = 1;

% LGN mosaic stuff from Greg's code 
filename = 'K042809001.nex'; 
CELLTYPE = 'P'; % example parvocell

% ecc_to_diam_deg = @(rf_r_deg)(0.0032*rf_r_deg+0.02)/sqrt(2); % Croner and Kaplan
HUMAN2MONKPSCALEFACTOR = .80; % From Dacey and Petersen. reasonable range: [.77 .81];
a = 0.9729; % Table 1
r2 = 1.084; % Table 1
re = 7.633; % Table 1
dc_0 = 14804.6; % Cone density of fovea
rm = 41.03; % See Equation 7
ecc_to_diam_deg = @(x)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
    (2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)))...
    ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
    *HUMAN2MONKPSCALEFACTOR); % Monkey midget RFs are slightly smaller than human midget RFs


% Support function for calculation of mus and S2
bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2));

stro = nex2stro(findfile(filename));
sigma_gabor = 0.15; % DVA
sigmas_n = 1;
rf_ecc = sqrt((stro.sum.exptParams.rf_x/10).^2+(stro.sum.exptParams.rf_y/10).^2);
RF_diam_deg = ecc_to_diam_deg(rf_ecc); % 2 SDs of the Gaussian RF
RFdistance = RF_diam_deg;
[x_deg,y_deg] = meshgrid(linspace(-1.0,1.0,11));
RFTRUNCATIONINSD = 1;

Rad3Over2 = sqrt(3)/2;
x_centers = [x_deg(1,1):RFdistance:x_deg(1,end)];
closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
x_centers = x_centers-x_centers(closest_to_zero);
x_centers_mat = repmat(x_centers,size(x_centers,2),1);
y_centers = x_centers;
x_centers_mat = x_centers_mat*Rad3Over2;
y_centers_mat = repmat(y_centers',1,size(y_centers,2));
y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RFdistance;
mus = zeros(numel(x_centers_mat),1);
interRFdistances = zeros(numel(x_centers_mat),numel(x_centers_mat));

% Loading the STA and the Gabor fits of the example file
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
load Gaborparams.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
% Singleopponent = logical(Singleopponent);
Singleopponent = logical(Ratio_of_power<1.2);

II = 132; % Index of an example cell
lambda = 12;

%STA
STA = Output_List{II,2};
normfactor = 0.5/(max(abs(STA(:)))+0.01);
STA = normfactor*STA + 0.5;
STA = reshape(STA,[10 10 3]);
STA = sigmoid(STA,11,0.5);
    
% SVD derived spatial RF 
imRF = Output_List{II,4};
imRF = sigmoid(imRF,lambda,0)-0.5;

% Gabor fit 
fit = modelfits{II,2};
fit = sigmoid(fit,lambda,0)-0.5;

% fft of the Gabor fit 
fitfft = fftshift(fft2(fit));

% Plotting the Gabor function
Normalizationfactor = 5;
params = Gaborparams{II};
theta = params.theta;
lambda = params.lambda/Normalizationfactor;
phi = params.phi;
sigma = params.sigma/Normalizationfactor;
gamma = params.gamma;
xoffset = params.xoffset/Normalizationfactor;
yoffset = params.yoffset/Normalizationfactor;
amplitude = params.amplitude;  

interval = [1:1:size(fit,1)]/Normalizationfactor;
interval = interval-median(interval);
[X, Y] = meshgrid(interval);
X = X-xoffset; Y = Y+yoffset; % So negative numbers mean down
xprime = X.*cos(-theta)+Y.*sin(-theta);
yprime = -X.*sin(-theta)+Y.*cos(-theta);
gabor = amplitude*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos((2.*pi.*yprime./lambda)-phi);
f = @(x,y) amplitude*exp(-((-1*(x-xoffset)*sin(-theta)+(y-yoffset)*cos(-theta))^2+gamma^2*((x-xoffset)*cos(-theta)+(y-yoffset)*sin(-theta))^2)/(2*sigma^2))*cos((2*pi*((x-xoffset)*cos(-theta)+(y-yoffset)*sin(-theta))/lambda)-phi);

%
figure(plot_counter); subplot(221); hold on;
tmp = (RF_diam_deg/2)*[cos(linspace(0,2*pi,200))', sin(linspace(0,2*pi,200))'];
for j = 1:numel(x_centers_mat)
    plot(x_centers_mat(j)+tmp(:,1),y_centers_mat(j)+tmp(:,2),'k-');
end
fcontour(f,'Linewidth',2); colormap(jet); axis square;
set(gca,'Tickdir','out','Xlim',[x_deg(1), x_deg(end)],'Ylim',[y_deg(1), y_deg(end)]); hold off;

figure(plot_counter); subplot(222); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; title('STA')
subplot(223); imagesc(255*(fit./(2*max(abs(fit(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]); title('Gabor fit');
subplot(224); imagesc(255*(imRF./(2*max(abs(imRF(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]); title('Spatial weighting function');

% Another plot for understanding the number of midget cells providing input to individual V1 cells
ECC = 0.5:0.5:7;
MIDGET_NUM = (1./ecc_to_diam_deg(ECC)).^2;

figure(plot_counter); plot(ECC,MIDGET_NUM,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out','XTick',1:1:7,'XScale','log'); xlabel('Eccentricity (deg)'); ylabel('Midget RGCs density');
plot_counter = plot_counter + 1;

%% Section 2: Simulating resulting color tuning with random L- and M-weighted input

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
% Singleopponent = logical(Singleopponent);
Singleopponent = logical(Ratio_of_power<1.2);

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];


% Implementing a bootstrap
TRIALS = 50000;
L_wt = []; 
M_wt = [];
S_wt = [];
for ii = 1:TRIALS
    S = 0.7;
%     while S>0.2
        tmp = randn(1,3);
        tmp = sign(tmp(2))*tmp./sum(abs(tmp));
        S = abs(tmp(3));
%     end
    L_wt = [L_wt; tmp(1)];
    M_wt = [M_wt; tmp(2)];
    S_wt = [S_wt; S];
    
end

% Plotting the cone weights 
figure(plot_counter); set(gcf,'Name','Cone wts') 
subplot(121); plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');

subplot(122); plot(L_wt,M_wt,'k.'); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); 
plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M'); hold off;
plot_counter = plot_counter + 1;

% Comparing the bootstrapped and the true distribution
figure(plot_counter); 

subplot(221); histogram(conewts_svd(3,[LumIds_conewts])+conewts_svd(2,[LumIds_conewts]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(conewts_svd(3,[ColorOpponentIds_conewts])+conewts_svd(2,[ColorOpponentIds_conewts]),-1:0.1:1,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(conewts_svd(3,[Other_conewts])+conewts_svd(2,[Other_conewts]),-1:0.1:1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
histogram(conewts_svd(3,[Sconedominated_conewts])+conewts_svd(2,[Sconedominated_conewts]),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0.5 1],'Linewidth',2);
axis square; set(gca,'Tickdir','out'); xlabel('S+M signal'); ylabel('P');

subplot(222); histogram(conewts_svd(1,[LumIds_conewts])+conewts_svd(2,[LumIds_conewts]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(conewts_svd(1,[ColorOpponentIds_conewts])+conewts_svd(2,[ColorOpponentIds_conewts]),-1:0.1:1,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(conewts_svd(1,[Other_conewts])+conewts_svd(2,[Other_conewts]),-1:0.1:1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
histogram(conewts_svd(1,[Sconedominated_conewts])+conewts_svd(2,[Sconedominated_conewts]),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0.5 1],'Linewidth',2);
axis square; set(gca,'Tickdir','out'); xlabel('L+M signal'); ylabel('P');

subplot(223); histogram(conewts_svd(3,[LumIds_conewts ColorOpponentIds_conewts Other_conewts Sconedominated_conewts])+conewts_svd(2,[LumIds_conewts ColorOpponentIds_conewts Other_conewts Sconedominated_conewts]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'Normalization','probability'); hold on;
histogram(S_wt+M_wt,-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2,'Normalization','probability'); hold on;
axis square; set(gca,'Tickdir','out'); xlabel('S+M signal'); ylabel('P');

subplot(224); histogram(conewts_svd(1,[LumIds_conewts ColorOpponentIds_conewts Other_conewts Sconedominated_conewts])+conewts_svd(2,[LumIds_conewts ColorOpponentIds_conewts Other_conewts Sconedominated_conewts]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'Normalization','probability'); hold on;
histogram(L_wt+M_wt,-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2,'Normalization','probability'); hold on;
axis square; set(gca,'Tickdir','out'); xlabel('L+M signal'); ylabel('P');

plot_counter = plot_counter + 1;


