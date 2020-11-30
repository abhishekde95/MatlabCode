%% This script is for generating Fig 1 for a Vision Research paper

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
gl.nRnds = 3;
gl.nPres = 1;
gl.colmaxcc = .09;
gl.lumscale = 4;
gl.ang = 0;
params0 = [50 gl.colmaxcc 100 3 0 gl.ang];

% Construct Polar Grid
if mod(gl.nRnds,2) == 0
    thetaspace = thetaspace / 2^((gl.nRnds-2)/2);
    rhospace = rhospace / 2^(gl.nRnds/2);
elseif mod(gl.nRnds,2) == 1
    thetaspace = thetaspace / 2^((gl.nRnds-1)/2);
    rhospace = rhospace / 2^((gl.nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
Lcc = tempLcc .* gl.colmaxcc;
Mcc = tempMcc .* gl.colmaxcc;
stim_RS = [Lcc Mcc];

% Scale points linearly
linscalemat = [gl.lumscale 0; 0 1] * [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
stim_LinDist = stim_RS * linscalemat;

% Whitened points (RS points + rotation to distinguish from RS)
arbrot = 0;
whitemat = inv(linscalemat) * [cos(arbrot) sin(arbrot); -sin(arbrot) cos(arbrot)];
stim_LinDistWhite = stim_LinDist * whitemat;

% Nonlinear scaling of stimuli (keeps angles constant)
gl.lummaxcc = gl.colmaxcc*gl.lumscale;
nlscale = gl.lummaxcc*gl.colmaxcc./sqrt((gl.colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(gl.lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* nlscale;
Mcc = tempMcc .* nlscale;
stim_NLDist = [Lcc Mcc];



%% Spike-Triggered Average Stimulus

% Responses to stimuli
nsp_RS = ComputeNakaRushtonJPW(params0,stim_RS,'surface7');
nsp_LinDist = ComputeNakaRushtonJPW(params0,stim_LinDist,'surface7');
nsp_NLDist = ComputeNakaRushtonJPW(params0,stim_NLDist,'surface7');

% Calculate STA in all spaces
STA_RS = (nsp_RS./sum(nsp_RS))' * stim_RS;
STA_LinDist = (nsp_LinDist./sum(nsp_LinDist))' * stim_LinDist;
STA_LinDistWhite = (nsp_LinDist./sum(nsp_LinDist))' * stim_LinDistWhite;
STA_NLDist = (nsp_NLDist./sum(nsp_NLDist))' * stim_NLDist;


%% Set Up Some Figure Stuff

% Surface for RS space
x = linspace(min(stim_RS(:,1))*1.1,max(stim_RS(:,1))*1.1,100);
y = linspace(min(stim_RS(:,2))*1.1,max(stim_RS(:,2))*1.1,100);
[X_RS,Y_RS] = meshgrid(x,y);
surf = ComputeNakaRushtonJPW(params0,[X_RS(:) Y_RS(:)],'surface7');
surf_RS = reshape(surf,size(X_RS));

% Surface for LinDist space
x = linspace(min(stim_LinDist(:,1))*1.1,max(stim_LinDist(:,1))*1.1,100);
y = linspace(min(stim_LinDist(:,2))*1.1,max(stim_LinDist(:,2))*1.1,100);
[X_LinDist,Y_LinDist] = meshgrid(x,y);
surf = ComputeNakaRushtonJPW(params0,[X_LinDist(:) Y_LinDist(:)],'surface7');
surf_LinDist = reshape(surf,size(X_LinDist));

% Surface for Whitened LinDist space
x = linspace(min(stim_LinDistWhite(:,1))*1.1,max(stim_LinDistWhite(:,1))*1.1,100);
y = linspace(min(stim_LinDistWhite(:,2))*1.1,max(stim_LinDistWhite(:,2))*1.1,100);
[X_LinDistWhite,Y_LinDistWhite] = meshgrid(x,y);
temppts = [X_LinDistWhite(:),Y_LinDistWhite(:)] * inv(whitemat);
surf = ComputeNakaRushtonJPW(params0,temppts,'surface7'); % Using responses from LinDist
surf_LinDistWhite = reshape(surf,size(X_LinDistWhite));

% Surface for NLDist
x = linspace(min(stim_NLDist(:,1))*1.1,max(stim_NLDist(:,1))*1.1,100);
y = linspace(min(stim_NLDist(:,2))*1.1,max(stim_NLDist(:,2))*1.1,100);
[X_NLDist,Y_NLDist] = meshgrid(x,y);
surf = ComputeNakaRushtonJPW(params0,[X_NLDist(:) Y_NLDist(:)],'surface7');
surf_NLDist = reshape(surf,size(X_NLDist));

%% Bubble Plot with Contours

% RS Distribution
figure(1); clf; hold on; box on;
set(gcf,'Name','Rotationally Symmetric Distribution')
bubblesize = (nsp_RS-min(nsp_LinDist))/(max(nsp_LinDist)-min(nsp_LinDist))*25+5;
contlevs = linspace(0,max(params0(1)),9);
for i = 1:size(stim_RS,1)
    h = plot3(stim_RS(i,1),stim_RS(i,2),nsp_RS(i),'ko'); hold on;
    set(h,'MarkerFaceColor','black','MarkerSize',bubblesize(i),'MarkerEdgeColor','white');
end
contour(X_RS,Y_RS,surf_RS,contlevs);
plot3(STA_RS(1),STA_RS(2),mean(nsp_RS),'co','MarkerFaceColor','c','MarkerSize',max(bubblesize)); % STA
[x,y] = pol2cart(gl.ang,1);
plot([0 x],[0 y],'r'); % true direction
lim = max(stim_RS(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);


%Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2a','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2a','-depsc');
% end
% disp('Fig 2a done.')



% LinDist Distribution
figure(2); clf; hold on; box on; axis equal;
set(gcf,'Name','Linearly Distended Distribution')
bubblesize = (nsp_LinDist-min(nsp_LinDist))/(max(nsp_LinDist)-min(nsp_LinDist))*25+5;
for i = 1:size(stim_LinDist,1)
    h = plot(stim_LinDist(i,1),stim_LinDist(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',bubblesize(i),'MarkerEdgeColor','white');
end
contour(X_LinDist,Y_LinDist,surf_LinDist,contlevs);
plot(STA_LinDist(1),STA_LinDist(2),'co','MarkerFaceColor','c','MarkerSize',max(bubblesize)); % STA
[x,y] = pol2cart(gl.ang,1);
plot([0 x],[0 y],'r'); % true direction
lim = max(stim_LinDist(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);

% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig2b','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig2b','-depsc');
% end
% disp('Fig 2b done.')

% Whitened LinDist Distribution
figure(3); clf; hold on; box on; axis equal;
set(gcf,'Name','Whitened Lineraly Distended Distribution')
for i = 1:size(stim_LinDistWhite,1)
    h = plot(stim_LinDistWhite(i,1),stim_LinDistWhite(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',bubblesize(i),'MarkerEdgeColor','white'); %bubblesize preserved from LinDist
end
contour(X_LinDistWhite,Y_LinDistWhite,surf_LinDistWhite,contlevs);
plot(STA_LinDistWhite(1),STA_LinDistWhite(2),'co','MarkerFaceColor','c','MarkerSize',max(bubblesize)); % STA
truedir = [x y] * inv(whitemat');
plot([0 truedir(1)],[0 truedir(2)],'r'); % true direction
lim = max(stim_LinDistWhite(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);

% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig3c','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig3c','-depsc');
% end
% disp('Fig 3c done.')


% NLDist Distribution
figure(4); clf; hold on; box on; axis equal;
set(gcf,'Name','Nonlinearly Distended Distribution')
bubblesize = (nsp_NLDist-min(nsp_LinDist))/(max(nsp_LinDist)-min(nsp_LinDist))*25+5;
for i = 1:size(stim_NLDist,1)
    h = plot(stim_NLDist(i,1),stim_NLDist(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',bubblesize(i),'MarkerEdgeColor','white');
end
contour(X_NLDist,Y_NLDist,surf_NLDist,contlevs);
plot(STA_NLDist(1),STA_NLDist(2),'co','MarkerFaceColor','c','MarkerSize',max(bubblesize)); % STA
[x,y] = pol2cart(gl.ang,1);
plot([0 x],[0 y],'r'); % true direction
lim = max(stim_NLDist(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);

disp('Not Saving figure 4.')