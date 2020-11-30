% For plotting and formatting data for Vision Reasearch Paper Figure 3.
% Modeled data are calculated in seperate scripts.

%% Figure 2a: RS Dataset Fitting Error

% Load variables
if ismac
    load('/Users/jpatrickweller/Dropbox/VisionResearchPaper//Fig2aData(RS).mat')
else
    load('C:\Users\jpweller\Dropbox\VisionResearchPaper\/Fig2aData(RS).mat')
end
fig2a = data;

% Set up figure
figure(1); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.1 .3 .8 .4])
plot([min(fig2a.angs) max(fig2a.angs)],[0 0],'k')
xlim([min(fig2a.angs) max(fig2a.angs)])
ylim([-max(abs(fig2a.RWA.means+fig2a.RWA.std))*1.1 max(abs(fig2a.RWA.means+fig2a.RWA.std))*1.1])

% Plot RWA data
h1 = shadedErrorBar(fig2a.angs,fig2a.RWA.means,fig2a.RWA.std);
set(h1.patch,'FaceColor',[0.5 0 0.5],'FaceAlpha',.1)

% Plot LL data
h2 = shadedErrorBar(fig2a.angs,fig2a.LL.means,fig2a.LL.std,'g');
set(h2.patch,'FaceAlpha',.5)

% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2a','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2a','-depsc');
% end
% disp('Fig 2a done.')


%% Figure 2b (Linearly Transformed RS Dataset)

% Load variables
if ismac
    load('/Users/jpatrickweller/Dropbox/VisionResearchPaper/FigbData(transRS).mat')
else
    load('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2bData(transRS).mat')
end
fig2b = data;

% Set up figure
figure(2); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.1 .3 .8 .4])
plot([min(fig2b.angs) max(fig2b.angs)],[0 0],'k')
xlim([min(fig2b.angs) max(fig2b.angs)])
maxx = max(abs(fig2b.RWAraw.means+fig2b.RWAraw.std))*1.1;
%ylim([-maxx maxx])
ylim([-65 65])

% Plot raw RWA data
h1 = shadedErrorBar(fig2b.angs,fig2b.RWAraw.means,fig2b.RWAraw.std,'r');
set(h1.patch,'FaceAlpha',.5)

% Plot corrected RWA data
h2 = shadedErrorBar(fig2b.angs,fig2b.RWAcorrect.means,fig2b.RWAcorrect.std,'b');
set(h2.patch,'FaceAlpha',.5)

% Plot LL data
h3 = shadedErrorBar(fig2b.angs,fig2b.LL.means,fig2b.LL.std,'g');
set(h3.patch,'FaceAlpha',.5)

% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2b','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2b','-depsc');
% end
% disp('Fig 2b done.')

%% Figure 2c (Nonlinearly Transformed RS Dataset)

% Load variables
if ismac
    load('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2cData(NLtransRS).mat')
else
    load('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2cData(NLtransRS).mat')
end
fig2c = data;

% Set up figure
figure(3); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.1 .3 .8 .4])
plot([min(fig2c.angs) max(fig2c.angs)],[0 0],'k')
xlim([min(fig2c.angs) max(fig2c.angs)])
maxx = max(abs(fig2c.RWAraw.means+fig2c.RWAraw.std))*1.1;
%ylim([-maxx maxx])
ylim([-65 65])

% Plot raw RWA data
h1 = shadedErrorBar(fig2c.angs,fig2c.RWAraw.means,fig2c.RWAraw.std,'r');
set(h1.patch,'FaceAlpha',.5)

% Plot corrected RWA data
h2 = shadedErrorBar(fig2c.angs,fig2c.RWAcorrect.means,fig2c.RWAcorrect.std,'b');
set(h2.patch,'FaceAlpha',.5)

% Plot LL data
h3 = shadedErrorBar(fig2c.angs,fig2c.LL.means,fig2c.LL.std,'g');
set(h3.patch,'FaceAlpha',.5)

% Format and save fig
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2c','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2c','-depsc');
% end
% disp('Fig 2c done.')


%% Figure 2d: RS distribution

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
nRnds = 3;
nPres = 1;
maxcc = .09;

% Construct Polar Grid
if mod(nRnds,2) == 0
    thetaspace = thetaspace / 2^((nRnds-2)/2);
    rhospace = rhospace / 2^(nRnds/2);
elseif mod(nRnds,2) == 1
    thetaspace = thetaspace / 2^((nRnds-1)/2);
    rhospace = rhospace / 2^((nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
Lcc = tempLcc .* maxcc;
Mcc = tempMcc .* maxcc;

% Plot figure
figure(4); clf; hold on; axis equal; box on;
h = plot(Lcc,Mcc,'ko');
set(h,'MarkerFaceColor','k')
xlim([min(Lcc)*1.1 max(Lcc)*1.1])
ylim([min(Mcc)*1.1 max(Mcc)*1.1])

% Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2d','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2d','-depsc');
end
disp('Fig 2d done.')


%% Figure 2e: Linearly Transformed RS Distribution

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
nRnds = 3;
nPres = 1;
colmaxcc = .09;
linscale = 4;

% Construct Polar Grid
if mod(nRnds,2) == 0
    thetaspace = thetaspace / 2^((nRnds-2)/2);
    rhospace = rhospace / 2^(nRnds/2);
elseif mod(nRnds,2) == 1
    thetaspace = thetaspace / 2^((nRnds-1)/2);
    rhospace = rhospace / 2^((nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
Lcc = tempLcc .* colmaxcc;
Mcc = tempMcc .* colmaxcc;

% Linear transform of RS Dataset
scalemat = [linscale 0; 0 1] * [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
newpts = [Lcc Mcc] * scalemat;
Lcc = newpts(:,1);
Mcc = newpts(:,2);

% Plot figure
figure(5); clf; hold on; axis equal; box on;
h = plot(Lcc,Mcc,'ko');
set(h,'MarkerFaceColor','k')
xlim([min(Lcc)*1.1 max(Lcc)*1.1])
ylim([min(Mcc)*1.1 max(Mcc)*1.1])

% Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2e','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2e','-depsc');
end
disp('Fig 2e done.')



%% Figure 2f: Nonlinearly Transformed RS Distribution

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
nRnds = 3;
nPres = 1;
colmaxcc = .09;
NLscale = 4;

% Construct Polar Grid
if mod(nRnds,2) == 0
    thetaspace = thetaspace / 2^((nRnds-2)/2);
    rhospace = rhospace / 2^(nRnds/2);
elseif mod(nRnds,2) == 1
    thetaspace = thetaspace / 2^((nRnds-1)/2);
    rhospace = rhospace / 2^((nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Linear transform of RS Dataset
[x,~] = pol2cart(3*pi/4,colmaxcc);
lumxy = abs(x * NLscale);
[~,lummaxcc] = cart2pol(lumxy,lumxy);
scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(lummaxcc.*sin(thetas-pi/4)).^2);
newpts = [tempLcc tempMcc] .* repmat(scale,[1 2]);
Lcc = newpts(:,1);
Mcc = newpts(:,2);

% Plot figure
figure(6); clf; hold on; axis equal; box on;
h = plot(Lcc,Mcc,'ko');
set(h,'MarkerFaceColor','k')
xlim([min(Lcc)*1.1 max(Lcc)*1.1])
ylim([min(Mcc)*1.1 max(Mcc)*1.1])

% Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2f','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2f','-depsc');
end
disp('Fig 2f done.')

