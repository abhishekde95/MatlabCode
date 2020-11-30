% Figures for thesis

% For making GLMS figures
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
datatypes = GLMSPopData(1,:);



%% Making 2D surface

idx = 135; % pan color ***
%idx = 44; % red horseshoe **
%idx = 63; % green horsehoe


GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};

%GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
%params = poppanel.params(idx,:);
params = poppanel.oneD.bichrom.parvals(idx,:);
params1 = params;
params1(4) = 0;
params1(3) = params1(2);
params2 = params1;
params2(2:3) = params(4);
params2(end-1) = params(end-1)+pi/2;

params = [100 1/.3 0 1/.5 3 1 pi 1];
%params1 = [100 1/.3 1/0 0 3 1 0 1];
%params2 = [100 1/.5 1/0 0 3 1 pi/2 1];

% Bubbble plot surface (this is annoying, but bubbles, contours, and
% surface must all be saved sepeartely for Illustrator to recognize it.
%meannsp = GLMP.subunit{sub}.meannspikes;
%ticks = -.08:.04:.08; % RS axes
ticks = -1:.2:1; % regular axes
cmap = repmat(linspace(.7,0,32)',1,3);

% Define thetas and rhos for the whole grid
% thetas = linspace(-pi,pi,361)';
% majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
% minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*2;
% nom = majorax * minorax;
% denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
% rhos = nom ./ sqrt(denom);
% scalars = 0:.01:1;
% rhosgrid = rhos * scalars;
% thetasgrid = repmat(thetas,1,numel(scalars));
%[x,y] = pol2cart(thetasgrid,rhosgrid);

xx = linspace(-1,1,101);
[x,y] = meshgrid(xx,xx);


surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],poppanel.surftype);
surface = reshape(surface,size(x));
% surface1 = ComputeNakaRushtonJPW(params1,[x(:) y(:)],poppanel.surftype);
% surface1 = reshape(surface1,size(x));
% surface2 = ComputeNakaRushtonJPW(params2,[x(:) y(:)],poppanel.surftype);
% surface2 = reshape(surface2,size(x));
axlim = max(x(:));

%%% (1) contours of 1 linear filter %%%
% figure(31); clf;  colormap(cmap);
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'renderer','paint')
% contour(x,y,surface1,5,'linewidth',2); hold on; % masked grid
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% name = 'contours1';
% box on; grid off;
% axis equal square tight
% set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%%% (2) contours of second linear filter %%%
% figure(32); clf;  colormap(cmap);
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'renderer','paint')
% contour(x,y,surface2,5,'linewidth',2); hold on; % masked grid
% %plot([0 PD(1)],[0 PD(2)],'k-')
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% name = 'contours2';
% box on; grid off;
% axis equal square tight
% set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%%% (2) contours of second linear filter %%%
% figure(32); clf;  colormap(cmap);
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'renderer','paint')
% contour(x,y,surface1+surface2,5,'linewidth',2); hold on; % masked grid
% %plot([0 PD(1)],[0 PD(2)],'k-')
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% name = 'contours3';
% box on; grid off;
% axis equal square tight
% set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% %print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


% %%% (3) Surface %%%
% figure(30); clf; colormap(cmap) 
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'renderer','paint')
% h = surf(x,y,surface); hold on;
% set(h,'edgecolor','none');
% alpha(.3);
% name = 'surf';
% box on; grid off;
% axis equal square tight
% set(gca,'XTick',[],'YTick',[],'CameraPosition',[0 0 275.3526],...
%     'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%%% (4) contours and pd %%%
figure(33); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
%plot([0 PD(1)],[0 PD(2)],'k-')
%set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Making just -lum alone bubbleplot
% idx = 38;
idx = 77; % **

bubblemaker(idx);

%% normalized log likelihood line for example neurons

idx1d = 70; % green unichrom **
idx2d = 135; % pan color **


ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
%params = poppanel.params(idx,:);


b1dgr = poppanel.oneD.bichrom.normLLs(idx1d);
u1dgr = poppanel.oneD.unichrom.normLLs(idx1d);
b2dgr = poppanel.twoD.bichrom.normLLs(idx1d);
u2dgr = poppanel.twoD.unichrom.normLLs(idx1d);

b1dpc = poppanel.oneD.bichrom.normLLs(idx2d);
u1dpc = poppanel.oneD.unichrom.normLLs(idx2d);
b2dpc = poppanel.twoD.bichrom.normLLs(idx2d);
u2dpc = poppanel.twoD.unichrom.normLLs(idx2d);

figure(1); clf; hold on;
set(gcf,'renderer','paint');
xlim([-.5 1.5])
ylim([-1.5 1.5])
plot([0 1],[0 0],'k-')
plot([0 0],[-1 1],'k-')
plot([1 1],[-1 1],'k-')

plot([b1dgr b1dgr],[-.5 .5],'b-')
plot([b2dgr b2dgr],[-.5 .5],'b-')

plot([b1dpc b1dpc],[-.5 .5],'r-')
plot([b2dpc b2dpc],[-.5 .5],'r-')

name = 'normll';
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')



%% 1D surfaces (unipolar and bipolar)

% params
A = 50;
bl = 0;
exp = 3;
sig1 = 1/.5;
sig2 = 0;
%sig2 = sig1*.65;
orthsig = 0;
rot = -pi/4;
params = [A sig1 sig1 orthsig exp bl rot];
cmap = repmat(linspace(.7,0,32)',1,3);


% Define thetas and rhos for the whole grid
thetas = linspace(-pi,pi,361)';
majorax = .9*1.1;
minorax = .09*1.1;
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


% surace
% x = -1:.01:1;
% [xx,yy] = meshgrid(x,x);
% surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
% surface = reshape(surface,size(xx));
% %surface = xx+yy;
% axlim = max(xx(:));

% plot Surface
figure(30); clf; colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
h = surf(x,x,surface); hold on;
set(h,'edgecolor','none');
alpha(.3);
name = 'surf';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'XTick',[],'YTick',[],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% plot contour
figure(31); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
contour(x,x,surface,5,'linewidth',2); hold on; % masked grid
name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'XTick',[],'YTick',[],'CameraPosition',[0 0 275.3526],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')




%% Figure 5A %%
% 1D polar PD histogram, smoothed polar hist, color modes, color wedges

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

%L = poppanel.oneDL;
%L = poppanel.oneDL & poppanel.monk.NutL;
L = poppanel.oneDL & poppanel.monk.MauiL;
angs = poppanel.tuning(L);

% Create histogram and smooth
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
PSTH = hist(angs,bins);
histangs = linspace(-pi,pi,21);

figure(31); clf; 
h = polarhistogram(angs,histangs,'Facecolor',[.5 .5 .5],'edgecolor','k'); hold on;
histmax = max(h.Values);
%h.DisplayStyle = 'stairs';
rlim([0 histmax*1.5])

name = 'mauipdhist';
%name = 'nutpdhist';
set(31,'PaperPositionMode','auto','renderer','paint')
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Figure 7A %%
% 2D polar PD histogram, smoothed polar hist, color modes, color wedges

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

ndL = poppanel.twoDL;
pcL = ndL & poppanel.twoD.bichrom.eliL;
sampn = 1000;

angmat = nan(sampn,4);
sigmat = nan(sampn,4);

% preferred directions
%PDangs = poppanel.tuning(pcL);
%PDsigs = 1./poppanel.params(pcL,2);
PDangs = randi(360,sampn,1);
%PDsigs = ones(size(PDangs));
angmat(:,1) = PDangs;
%sigmat(:,1) = PDsigs;

% APD
%apdsigs = 1./poppanel.params(pcL,3);
apdangs = PDangs+pi;
if any(apdangs>pi)
    apdangs(apdangs>pi) = apdangs(apdangs>pi)-(2*pi);
end
angmat(:,2) = apdangs;
%sigmat(:,2) = apdsigs;
%h = plot(apdangs,apdsigs,'ko');

% Orthogonal
%osigs = 1./poppanel.params(pcL,4);
oangs1 = PDangs+(pi/2);
if any(oangs1>pi)
    oangs1(oangs1>pi) = oangs1(oangs1>pi)-(2*pi);
end
angmat(:,3) = oangs1;
%sigmat(:,3) = osigs;
%h = plot(oangs1,osigs,'ko');

oangs2 = oangs1+pi;
if any(oangs2 > pi)
    oangs2(oangs2>pi) = oangs2(oangs2>pi)-(2*pi);
end
angmat(:,4) = oangs2;
%sigmat(:,4) = osigs';
%h = plot(oangs2,osigs,'ko');

% arrange lowest to highest 
[sortangmat,idx] = sort(angmat,2);
%sortsigmat = nan(size(sortangmat));
% for n = 1:size(sortsigmat,1)
%     sortsigmat(n,:) = sigmat(n,idx(n,:));
% end
%angs = poppanel.tuning(ndL & ~poppanel.excludeL);

% Create histogram and smooth
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 45*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH = hist(sortangmat(:),bins);
smoothPSTH = conv(PSTH,gaussfilter,'same');
smoothPSTH = smoothPSTH./max(smoothPSTH) .* max(PSTH); % normalizing by max height of PSTH for plotting purposes
histangs = linspace(-pi,pi,21);

% Plot histogram  
figure(60); clf;
h = polarhistogram(sortangmat,histangs,'FaceColor',[.5 .5 .5],...
    'edgecolor','k','linewidth',2); hold on;
histmax = max(h.Values);
rlim([0 histmax*1.5]);
thetalim([0 90]);

name = 'Figure 7A';
set(60,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Figure 9B %%%
% 2D PD and c50 means and standard deviations

%GLMSPopGUI_Params

cols = {'+lum' '-lum' 'red chrom' 'green chrom'};
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.twoDL;
radpm = 22.5/180*pi;
name = 'Figure 9B';

% Pull out pan color cells
hsL = ndL & poppanel.twoD.unichrom.eliL;
prinaxmat = nan(sum(hsL),4);
prinsigmat = nan(sum(hsL),4);

% preferred directions
PDangs = poppanel.tuning(hsL)./pi*180;
PDsigs = 1./poppanel.params(hsL,2);
prinaxmat(:,1) = PDangs;
prinsigmat(:,1) = PDsigs;

% APD
apdsigs = 1./poppanel.params(hsL,3);
apdangs = PDangs+180;
if any(apdangs>180)
    apdangs(apdangs>180) = apdangs(apdangs>180)-360;
end
prinaxmat(:,2) = apdangs;
prinsigmat(:,2) = apdsigs;

% Orthogonal
osigs = 1./poppanel.params(hsL,4);
oangs1 = PDangs+90;
if any(oangs1>180)
    oangs1(oangs1>180) = oangs1(oangs1>180)-360;
end
prinaxmat(:,3) = oangs1;
prinsigmat(:,3) = osigs;

oangs2 = PDangs-90;
if any(oangs2<-180)
    oangs2(oangs2<-180) = oangs2(oangs2<-180)+360;
end
prinaxmat(:,4) = oangs2;
prinsigmat(:,4) = osigs';

% organize from -pi to pi
[sortangmat,idx] = sort(prinaxmat,2);
sortsigmat = nan(size(sortangmat));
for n = 1:size(sortsigmat,1)
    sortsigmat(n,:) = prinsigmat(n,idx(n,:));
end

% Add padding to both sides (to emulate cyclical plotting)
paddedangmat = cat(2,sortangmat(:,1)-90,sortangmat,sortangmat(:,4)+90);
paddedsigmat = cat(2,sortsigmat(:,4),sortsigmat,sortsigmat(:,1));


% Figure 
figure(84); clf; hold on; box on;
set(gca,'xtick',linspace(-180,180,7),'tickdir','out')

% Color boxes, peaks inherrited from 8A
angL = nan(sum(hsL),4);
for n = 1:4
    xedge(1) = (peaks(n)+radpm)/pi*180;
    xedge(2) = (peaks(n)-radpm)/pi*180;
    angL(:,n) = PDangs < xedge(1) & PDangs > xedge(2);
    h = plot([xedge(1) xedge xedge(2) xedge(1)],[0 1 1 0 0]); % boxes
    set(h,'color',conpanel.carcols(n,:)); hold on;
end

% Calculate iso-c50 contours
linex = linspace(-pi,pi,181);
%linex = linspace(-pi,pi,9)
sigmat = nan(sum(hsL),numel(linex));
thetamat = nan(sum(hsL),numel(linex));
hsidx = find(hsL);
for n = 1:numel(hsidx)
    params = poppanel.params(hsidx(n),:);
    pdsig = 1/params(2);
    apdsig = 1/params(3);
    osig = 1/params(4);
    pd = params(end-1);
    shiftedx = linex-pd;
    shiftedx(shiftedx>pi) = shiftedx(shiftedx>pi)-(2*pi);
    shiftedx(shiftedx<-pi) = shiftedx(shiftedx<-pi) + (2*pi);
    posL = shiftedx<pi/2 & shiftedx>-pi/2;
    
    % Fit positive side with ellipse
    sigmat(n,posL) = (pdsig*osig) ./ sqrt((pdsig*sin(shiftedx(posL))).^2 + (osig*cos(shiftedx(posL))).^2);
    %sigmat(n,~posL) = abs(osig ./ sin(shiftedx(~posL)));
    sigmat(n,~posL) = nan;
    
end

% Smooth means and stds
x1 = nanmean(sigmat);
x2 = nanstd(sigmat);
order = 3;
framelen = 31;
smoothmean = sgolayfilt(x1,order,framelen);
smoothstd = sgolayfilt(x2,order,framelen);

% Plot smoothed
h = plot(linex./pi*180,sigmat,'k-');
%set(h,'color',[.85 .85 .85])
shadedErrorBar(linex./pi*180,smoothmean,smoothstd,'b')
%shadedErrorBar(linex./pi*180,nanmean(sigmat),smoothstd,'b')
%plot(linex./pi*180,nanmean(sigmat),'r')
alpha(.4)


%idx = 63; % green horsehoe
idx = 44; % red horseshoe **
%idx = 158; % +lum horseshsoe
%idx = 68; % -lum horseshoe
idxL = find(hsL)==idx;
plot(linex./pi*180,sigmat(idxL,:),'r-')


% wedge points
% L = any(angL,2);
% plot(PDangs(L),PDsigs(L),'ko','markerfacecolor','k')
% oangs = PDangs+90;
% oangs(oangs>180) = oangs(oangs>180)-360;
% plot(oangs(L),osigs(L),'ko')
% oangs = PDangs-90;
% oangs(oangs<-180) = oangs(oangs<-180)+360;
% plot(oangs(L),osigs(L),'ko')

% ungrouped pts
% L = ~any(angL,2);
% plot(PDangs(L),PDsigs(L),'o','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
% oangs = PDangs+90;
% oangs(oangs>180) = oangs(oangs>180)-360;
% plot(oangs(L),osigs(L),'o','color',[.5 .5 .5])
% oangs = PDangs-90;
% oangs(oangs<-180) = oangs(oangs<-180)+360;
% plot(oangs(L),osigs(L),'o','color',[.5 .5 .5])

% Figure housekeeping
xlim([-180 180]); ylim([0 1])
title('Horseshoe cells: principle axes and c50s')
set(gcf,'PaperPositionMode','auto','numbertitle','off','name',name);
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% RF locations and sizes, divided by LN/LNLN.



glmps = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
rfparams = GLMSPopData(2:end,strcmp(datatypes,'RF Params'));

% by 1D vs 2D
figure(1); clf; hold on;
plot(0,0,'r-')
plot(0,0,'b-')
plot(0,0,'kx')
for n = 1:numel(glmps)
    GLMP = glmps{n};
    rf = rfparams{n};
    
    if poppanel.oneDL(n)
        h = circle([GLMP.rf_x GLMP.rf_y],sqrt(rf.sqdva/2/pi),1000,'r-');
    else
        h = circle([GLMP.rf_x GLMP.rf_y],sqrt(rf.sqdva/2/pi),1000,'b-');
    end
   
end
axis equal
%plot(0,0,'xk')
legend('LN','LNLN')
set(gcf,'PaperPositionMode','auto','renderer','paint');
print('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/rf_by_model','-depsc')


% by monkey
figure(2); clf; hold on;
plot(0,0,'g-')
plot(0,0,'y-')
plot(0,0,'kx')
for n = 1:numel(glmps)
    GLMP = glmps{n};
    rf = rfparams{n};
    
    if poppanel.monk.NutL(n)
        h = circle([GLMP.rf_x GLMP.rf_y],sqrt(rf.sqdva/2/pi),1000,'g-');
    else
        h = circle([GLMP.rf_x GLMP.rf_y],sqrt(rf.sqdva/2/pi),1000,'y-');
    end   
   
end
axis equal
legend('Nut','Maui')
set(gcf,'PaperPositionMode','auto','renderer','paint');
print('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/rf_by_monkey','-depsc')

%% Distribution of Latency

glmps = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
latencies = nan(size(glmps));

for n = 1:numel(glmps)
    GLMP = glmps{n};
    latencies(n) = GLMP.countingwin(1);
end
bins = linspace(0,.1,21);

figure(1); clf; hold on;
h = histogram(latencies(poppanel.oneDL),bins);
h.DisplayStyle = 'stairs';
h.EdgeColor = 'r';

h = histogram(latencies(poppanel.twoDL),bins);
h.DisplayStyle = 'stairs';
h.EdgeColor = 'b';

ticks = [0:.01:.1];
name = 'latency dist';
set(gca,'tickdir','out','XTick',ticks)
set(gcf,'PaperPositionMode','auto','renderer','paint');
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

mean(latencies)
std(latencies)

