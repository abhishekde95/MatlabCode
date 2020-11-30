function [] = bubblemaker(idx)
global GLMSPopData


datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};

ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
params = poppanel.params(idx,:);


% Bubbble plot surface (this is annoying, but bubbles, contours, and
% surface must all be saved sepeartely for Illustrator to recognize it.
meannsp = GLMP.subunit{sub}.meannspikes;
%ticks = -.08:.04:.08; % RS axes
ticks = -.8:.2:.8; % regular axes
cmap = repmat(linspace(.7,0,32)',1,3);

% Define thetas and rhos for the whole grid
thetas = linspace(-pi,pi,361)';
majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
nom = majorax * minorax;
denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
rhos = nom ./ sqrt(denom);
scalars = 0:.01:1;
rhosgrid = rhos * scalars;
thetasgrid = repmat(thetas,1,numel(scalars));
[x,y] = pol2cart(thetasgrid,rhosgrid);
surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],poppanel.surftype);
surface = reshape(surface,size(rhosgrid));
axlim = max(x(:));

%%% (1) Surface %%%
figure(30); clf; colormap(cmap) 
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')

h = surf(x,y,surface); hold on;
set(h,'edgecolor','none');
alpha(.3);

name = 'surf';
box on; grid off;
axis equal square tight
set(gca,'XTick',[],'YTick',[],'CameraPosition',[0 0 275.3526],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%%% (2) contours and pd %%%
figure(31); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')

rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
PD = rotmat * [0 .5]';
contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
plot([0 PD(1)],[0 PD(2)],'k-')
set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%%% (3) Bubbles %%%
figure(32); clf;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
maxnsp = max(max(GLMP.subunit{sub}.meannspikes));

polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
scalar = 30;
for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
    mn = meannsp(i)/maxnsp*scalar+scalar/10;
    h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko'); 
    set(h,'MarkerFaceColor','none','MarkerSize',mn)
end
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

% legend
L = GLMP.subunit{sub}.meannspikes == max(GLMP.subunit{sub}.meannspikes);
h = polar(GLMP.subunit{sub}.uniquetheta(L),GLMP.subunit{sub}.uniquerho(L),'ko');
temp = (GLMP.subunit{sub}.meannspikes(L)/maxnsp);
set(h,'MarkerFaceColor','none','MarkerSize',temp*scalar+scalar/10)
L = GLMP.subunit{sub}.meannspikes == min(GLMP.subunit{sub}.meannspikes);
k = polar(GLMP.subunit{sub}.uniquetheta(L),GLMP.subunit{sub}.uniquerho(L),'ko');
set(k,'MarkerFaceColor','none','MarkerSize',scalar/10)
frs = num2str(cat(1,round(max(GLMP.subunit{sub}.meanfr)),round(min(GLMP.subunit{sub}.meanfr))));
legend([h;k],frs,'location','northwest')

name = 'bubbles';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


end
