


if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{2,strcmp(datatypes,'GLMP')};


%%
%funds = GLMP.fundamentals ./ repmat(sum(GLMP.fundamentals),size(GLMP.fundamentals,1),1)
%M = GLMP.monspds' * funds;
M = GLMP.monspds' * GLMP.fundamentals;

% Construct stim dist in RGB space
thetas = shiftdim(0:pi/8:15*pi/8);
rhos = shiftdim(0:.25:1);
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = [thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))];
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);
[Rgun,Ggun] = pol2cart(thetas,rhos);
Bgun = zeros(size(Rgun));

% cone contrast space
cc = [Rgun Ggun Bgun] * M;
Lcc = cc(:,1);
Mcc = cc(:,2);

% Generate responses
params = [50 max(Lcc)/2 100 3 0 0];
nsp = ComputeNakaRushtonJPW(params,[Lcc Mcc],'surface7');

% Generate RGB surface
x = linspace(-1.1,1.1,100);
[X_RGB,Y_RGB] = meshgrid(x,x);
xylms = [X_RGB(:) Y_RGB(:) zeros(size(X_RGB(:)))] * M;
X_lms = xylms(:,1);
Y_lms = xylms(:,2);
surf_lms = ComputeNakaRushtonJPW(params,[X_lms Y_lms],'surface7');
surf_RGB = reshape(surf_lms,size(X_RGB));

% Generate LMS surface
x = linspace(min(Lcc)*1.1,max(Lcc)*1.1,100);
[X_LMS,Y_LMS] = meshgrid(x,x);
surf_LMS = ComputeNakaRushtonJPW(params,[X_LMS(:) Y_LMS(:)],'surface7');
surf_LMS = reshape(surf_LMS,size(X_LMS));

% Plot gun space
figure(1); clf; hold on; box on; axis equal
set(gcf,'Name','Gun Space')
bubblesize = (nsp-min(nsp))/(max(nsp)-min(nsp))*25+5;
contlevs = linspace(0,max(nsp),9);
contour(X_RGB,Y_RGB,surf_RGB,contlevs)
for i = 1:numel(Rgun)
    h = plot(Rgun(i),Ggun(i),'ko'); hold on;
    set(h,'MarkerSize',bubblesize(i),'MarkerEdgeColor','white','MarkerFaceColor','black');
end
prefdir = ([max(X_LMS(:)) 0 0] * M').*1000;
plot([0 prefdir(1)],[0 prefdir(2)],'r'); % true direction
lim = max(X_RGB(:));
xlim([-lim lim]); ylim([-lim lim]);
ticks = linspace(-1,1,5);
set(gca,'xtick',ticks,'xticklabel',ticks',...
    'ytick',ticks,'yticklabel',ticks)

% Format and save figure
set(gcf,'PaperPositionMode','auto')
name = 'Fig4b_RV';
if ismac
    print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
else ispc
    print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
end
disp(['Fig ' name ' done.'])


% Plot cone space
figure(2); clf; axis equal; hold on; box on;
set(gcf,'Name','Cone Contrast Space')
bubblesize = (nsp-min(nsp))/(max(nsp)-min(nsp))*25+5;
contlevs = linspace(0,max(nsp),9);
contour(X_LMS,Y_LMS,surf_LMS,contlevs)
for i = 1:numel(Lcc)
    h = plot(Lcc(i),Mcc(i),'ko'); hold on;
    set(h,'MarkerSize',bubblesize(i),'MarkerEdgeColor','white','MarkerFaceColor','black');
end
plot([0 max(X_LMS(:))],[0 0],'r'); % true direction
lim = max(X_LMS(:));
xlim([-lim lim]); ylim([-lim lim]);
ticks = linspace(-.04,.04,5);
set(gca,'xtick',ticks,'xticklabel',ticks',...
   'ytick',ticks,'yticklabel',ticks)

% Format and save figure
set(gcf,'PaperPositionMode','auto')
name = 'Fig3a_RV';
if ismac
    print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
else ispc
    print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
end
disp(['Fig ' name ' done.'])
