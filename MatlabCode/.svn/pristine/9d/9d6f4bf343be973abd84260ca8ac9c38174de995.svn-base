% For building figure 2 for VR_Rrevised


% Construct stim dist in RGB space
thetas = shiftdim(0:pi/8:15*pi/8);
rhos = shiftdim(0:.25:1);
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = [thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))];
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);
[Lcc_RS,Mcc_RS] = pol2cart(thetas,rhos);
Scc_RS = zeros(size(Lcc_RS));

% Convert from RS to distended distribution
M = [1 0 0; 0 0 0; 0 0 0];
cc = [Lcc_RS Mcc_RS Scc_RS] * M;
Lcc_D = cc(:,1);
Mcc_D = cc(:,2);

% Generate responses
params = [50 max(Lcc_RS)/2 100 3 0 0];
nsp_RS = ComputeNakaRushtonJPW(params,[Lcc_RS Mcc_RS],'surface7');
nsp_D = ComputeNakaRushtonJPW(params,[Lcc_D Mcc_D],'surface7');

% Generate surface
x = linspace(-1,1,100);
[X,Y] = meshgrid(x,x);
surf = ComputeNakaRushtonJPW(params,[X(:) Y(:)],'surface7');
surf = reshape(surf,size(X));

% Plot RS distribution
figure(1); clf; hold on; box on; axis equal
set(gcf,'Name','Fig 2B')
bubblesize = (nsp_RS-min(nsp_RS))/(max(nsp_RS)-min(nsp_RS))*25+5;
contlevs = linspace(0,max(nsp_RS),9);
contour(X,Y,surf,contlevs)
for i = 1:numel(Lcc)
    h = plot(Lcc_RS(i),Mcc_RS(i),'ko'); hold on;
    set(h,'MarkerSize',bubblesize(i),'MarkerFaceColor','black');
end
plot([0 1],[0 0],'r'); % true direction
lim = max(X(:));
xlim([-lim lim]); ylim([-lim lim]);
ticks = linspace(-1,1,5);
set(gca,'xtick',ticks,'xticklabel',ticks',...
    'ytick',ticks,'yticklabel',ticks)

%% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% name = 'Fig3b_RV';
% if ismac
%     print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
% else ispc
%     print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
% end
% disp(['Fig ' name ' done.'])


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
