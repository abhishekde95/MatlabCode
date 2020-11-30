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
GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');


%% Figure 1a
% This figure is just a snapshot of the white noise stimulus

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\nex files\';
end
DN = GLMSPopData{end,strcmp(datatypes,'DN')};
rawdata = nex2stro([char(library) DN.datafile '.nex']);

lumCC = sqrt((max(DN.lumCC)^2)/2);
colCC = sqrt((max(DN.colCC)^2)/2);
ccs = [lumCC lumCC 0; -lumCC -lumCC 0; colCC -colCC 0; -colCC colCC 0];
nstix = max(DN.NStixGrid);
pixperstix = 15;
randorder = randi(4,nstix^2,1);
bkgndrgb = rawdata.sum.exptParams.bkgndrgb';
bkgndlms =  bkgndrgb * (DN.M);
stimlms = ccs .* repmat(bkgndlms,size(ccs,1),1) + repmat(bkgndlms,size(ccs,1),1);
stimrgb = stimlms * inv(DN.M);

WNframe = nan(nstix*pixperstix,nstix*pixperstix,3);
temp = stimrgb(randorder,:);
temp = reshape(temp,[nstix nstix 3]);
for m = 1:nstix
    for p = 1:nstix
        t = repmat(temp(m,p,:),[pixperstix pixperstix 1]);
        WNframe((pixperstix*(m-1)+1):(pixperstix*m),(pixperstix*(p-1)+1):(pixperstix*p),:) = t;
    end
end

figure(10); clf; 
set(10,'units','normalized','pos',[.3 .3 .4 .4])
oneA = axes('parent',10,'units','normalized','pos',[.1 .1 .8 .8]);
image(WNframe,'parent',oneA);
set(gca,'xtick',[],'ytick',[]); axis square;
name = 'Figure 1a';
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%%% Figure 1B %%%
    
WNframes2 = nan(nstix*pixperstix,nstix*pixperstix,3,4);
for n = 1:4
    stimrgb2 = repmat(bkgndrgb,4,1);
    stimrgb2(n,:) = stimrgb(n,:);
    temp = stimrgb2(randorder,:);
    temp = reshape(temp,[nstix nstix 3]);
    for m = 1:nstix
        for p = 1:nstix
            t = repmat(temp(m,p,:),[pixperstix pixperstix 1]);
            WNframes2((pixperstix*(m-1)+1):(pixperstix*m),(pixperstix*(p-1)+1):(pixperstix*p),:,n) = t;
        end
    end
end

name = 'pLpM';
figure(11); clf; 
set(11,'units','normalized','pos',[.3 .3 .4 .4])
image(WNframes2(:,:,:,1));
set(gca,'xtick',[],'ytick',[]); axis square;
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

name = 'mLmM';
figure(12); clf;
set(12,'units','normalized','pos',[.3 .3 .4 .4])
image(WNframes2(:,:,:,2));
set(gca,'xtick',[],'ytick',[]); axis square;
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

name = 'pLmM';
figure(13); clf;
set(13,'units','normalized','pos',[.3 .3 .4 .4])
image(WNframes2(:,:,:,3));
set(gca,'xtick',[],'ytick',[]); axis square;
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

name = 'mLpM';
figure(14); clf;
set(14,'units','normalized','pos',[.3 .3 .4 .4])
image(WNframes2(:,:,:,4));
set(gca,'xtick',[],'ytick',[]); axis square;
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Figure 1C %%%
% This figure is an STA of a simple cell

%idx = 78; % 2 subunit simple cel
%idx = 153; % -lum unichrom
%idx = 82; % +lum unichrom but currently excluded
%idx = 33; % red unichrom
idx = 70; % green unichrom **
%idx = 92; % bichromatic (top choice)
%idx = 43; %bichromatic
%idx = 68; % -lum horseshoe
%idx = 63; % green horsehoe
%idx = 44; % red horseshoe **
%idx = 159; % +lum horseshoe
%idx = 154; % pan color

cellsel = ['cellselect(' num2str(idx) ',[])'];
subonoff = 'subonoff';
GLMSPopGUI_Tuning(cellsel,subonoff)

tunfig = get(150,'UserData');
cellpan = get(tunfig.cellpanel,'userdata');


figure(1); clf;
pLpM = copyobj(cellpan.pLpM.axes,gcf);
set(gca,'units','normalized','pos',[.1 .1 .8 .8]); 
axis square;
set(1,'PaperPositionMode','auto','renderer','paint')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' 'pLpM'],'-depsc')

figure(2); clf;
pLmM = copyobj(cellpan.pLmM.axes,gcf);
set(gca,'units','normalized','pos',[.1 .1 .8 .8]); 
axis square;
set(2,'PaperPositionMode','auto','renderer','paint')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' 'pLmM'],'-depsc')

figure(3); clf;
mLmM = copyobj(cellpan.mLmM.axes,gcf);
set(gca,'units','normalized','pos',[.1 .1 .8 .8]); 
axis square;
set(3,'PaperPositionMode','auto','renderer','paint')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' 'mLmM'],'-depsc')

figure(4); clf;
mLpM = copyobj(cellpan.mLpM.axes,gcf);
set(gca,'units','normalized','pos',[.1 .1 .8 .8]);
axis square;
set(4,'PaperPositionMode','auto','renderer','paint')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' 'mLpM'],'-depsc')

% set(pLmM,'units','normalized','pos',[.6 .1 .35 .35]); axis(pLmM,'square');
% set(mLpM,'units','normalized','pos',[.1 .1 .35 .35]); axis(mLpM,'square');
% set(mLmM,'units','normalized','pos',[.1 .6 .35 .35]); axis(mLmM,'square');
% set(15,'PaperPositionMode','auto')
% name = 'Figure 1c';
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Figure 1D %%%
% This figure is a depiction of points in the LM plane. It should highlight
% both WN contrasts and GLMP contrasts.


% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
nRnds = 3;
lummaxcc = .9;
colmaxcc = .09;

% Construct Polar Grid
if mod(nRnds,2) == 0
    thetaspace = thetaspace / 2^((nRnds-2)/2);
    rhospace = rhospace / 2^(nRnds/2);
else
    thetaspace = thetaspace / 2^((nRnds-1)/2);
    rhospace = rhospace / 2^((nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = [thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))];
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = lummaxcc * colmaxcc ./ sqrt((colmaxcc .* cos(thetas-pi/4)).^2 ...
    +(lummaxcc .* sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
Scc = zeros(size(Lcc));

% convert to display colors
DN = GLMSPopData{end,strcmp(datatypes,'DN')};
%rawdata = nex2stro([char(library) 'nex files\' DN.datafile '.nex']);
bkgndrgb = rawdata.sum.exptParams.bkgndrgb';
bkgndlms =  bkgndrgb * (DN.M);
stimlms = [Lcc Mcc Scc] .* repmat(bkgndlms,size(Lcc,1),1) + repmat(bkgndlms,size(Lcc,1),1);
stimrgb = stimlms * inv(DN.M);

% Plot
figure(16); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.2 .1 .6 .8]);
%set(gca,'tickdir','out'); axis square
set(gca,'xtick',[],'ytick',[]); axis square
[t,r] = cart2pol(Lcc, Mcc);
[~,order] = sortrows([r t]);
%plot(0,0,'ok','MarkerFaceColor',bkgndrgb,'MarkerSize',16);
for n = 1:numel(order)
    h = plot(Lcc(order(n)),Mcc(order(n)),'square','MarkerFaceColor',stimrgb(order(n),:),...
        'MarkerEdgeColor','w','MarkerSize',16);
end
xlim([min(Lcc)*1.1 max(Lcc)*1.1]);
ylim([min(Mcc)*1.1 max(Mcc)*1.1]);

set(16,'PaperPositionMode','auto','renderer','paint')
name = 'stimdist';
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-eps')
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')




%% Figure 2 %%

% generic NKR
x = -1:.01:1;
A = 50;
bl = 0;
exp = 3;
sig1 = .5;
sig2 = sig1/.65;
orthsig = 1/1.5;
params = [A A sig1 sig2 exp exp bl];
y = ComputeNakaRushtonJPW(params,x,'asymmetric');
y2 = y;
y2(1:find(x==0)) = 0;

figure(20); clf; hold on; box on;
name = 'nakarushton';
plot(x(1:find(x==0)),y(1:find(x==0)),'k--');
plot(x(1:find(x==0)),y2(1:find(x==0)),'k-');
plot(x(find(x==0):end),y(find(x==0):end),'k-')
xlim([min(x) max(x)])
ylim([-5 max(y)])
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-eps')

%%%  exponentials %%%
% regular old squaring normal style
x = -1:.01:1;
y = x.^2;
name = 'exp1';
figure(21); clf; hold on; box on;
plot(x(1:find(x==0)),y(1:find(x==0)),'k-');
plot(x(find(x==0):end),y(find(x==0):end),'k-')
xlim([min(x) max(x)])
ylim([-.1 max(y)])
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-eps')

% half squaring
x = -1:.01:1;
y = x.^2;
y2 = y;
y2(1:find(x==0)) = 0;
name = 'exp0';
figure(22); clf; hold on; box on;
plot(x(1:find(x==0)),y(1:find(x==0)),'k--');
plot(x(1:find(x==0)),y2(1:find(x==0)),'k--');
plot(x(find(x==0):end),y(find(x==0):end),'k-')
xlim([min(x) max(x)])
ylim([-.1 max(y)])
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-eps')

% reguar squaring
x = -1:.01:1;
y = x.^2;
name = 'exp2';
figure(23); clf; hold on; box on;
plot(x(1:find(x==0)),y(1:find(x==0)),'k--');
plot(x(find(x==0):end),y(find(x==0):end),'k-')
xlim([min(x) max(x)])
ylim([-.1 max(y)])
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
%export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-eps')

%% 1D surfaces (unipolar and bipolar)

% params
A = 50;
bl = 0;
exp = 3;
sig1 = 1/.5;
sig2 = 0;
%sig2 = sig1*.65;
orthsig = 0;
rot = pi/4;
params = [A sig1 sig2 orthsig exp bl rot];

% surace
x = -1:.01:1;
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
surface = reshape(surface,size(xx));
axlim = max(xx(:));
cmap = repmat(linspace(.7,0,32)',1,3);

% Surface
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
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% contour
figure(31); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')

rotmat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
PD = rotmat * [.5 0]';
contour(x,x,surface,5,'linewidth',2); hold on; % masked grid
%plot([0 PD(1)],[0 PD(2)],'k-')

name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'XTick',[],'YTick',[],'CameraPosition',[0 0 275.3526],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% 2D surface %%%

% Set params
A = 50;
bl = 1;
exp = 3;
sig1 = .5;
%sig2 = inf; osig = -.5; % unipolar hypertuned
%sig2 = sig1*1.65; osig = -.5; % bipolar hypertuned
sig2 = Inf; osig = .75; % horseshoe
%sig2 = sig1; osig = .75; % pan color
rot = -pi/4;
params = [A 1/sig1 1/sig2 1/osig exp bl rot];

% surface
x = -1:.01:1;
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
surface = reshape(surface,size(xx));
axlim = max(xx(:));
cmap = repmat(linspace(.7,0,32)',1,3);

% Surface
figure(30); clf; colormap(cmap)
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
h = surf(xx,yy,surface); hold on;
set(h,'edgecolor','none');
alpha(.3);

name = 'surf';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'XTick',[],'YTick',[],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% contour
figure(31); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')

rotmat = [cos(rot) -sin(rot); sin(rot) cos(rot)];
PD = rotmat * [.75 0]';
oax = rotmat * [0 .75; 0 -.75]';
contour(x,x,surface,5,'linewidth',2); hold on; % masked grid
plot([0 PD(1)],[0 PD(2)],'k-')
plot(oax(1,:),oax(2,:),'k-')

name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'XTick',[],'YTick',[],'CameraPosition',[0 0 275.3526],...
    'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')



%% bichromatic hyperbola %%%
name = 'hyperbola_bi';
osig = -2*orthsig;
params = [A sig1 sig2 osig exp bl rot];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz0 = reshape(z,size(xx));
params = [A sig1 0 osig exp bl rot];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz1 = reshape(z,size(xx));
params = [A sig2 0 osig exp bl rot+pi];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz2 = reshape(z,size(xx));

figure(25); clf; hold on; box on; axis equal;
set(gcf,'units','pixels','pos',[200 200 500 500])
surface(xx,yy,zz0,'edgealpha',0,'facealpha',.3)
contour3(xx,yy,zz1,'linewidth',2);
contour3(xx,yy,zz2,'linewidth',2,'linestyle','--');
colormap cool
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-jpg')

%% unichromatic hyperbola%%%
name = 'hyperbola_uni';
[xx,yy] = meshgrid(x,fliplr(x));
sig1 = 1/.5;
sig2 = 0;
osig = -2*orthsig;
bl = 0;
rot = -pi/4;
params = [A sig1 sig2 osig exp bl rot];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz0 = reshape(z,size(xx));

figure(24); clf; hold on; box on; axis equal;
set(gcf,'units','pixels','pos',[200 200 500 500])
surf(xx,yy,zz0,'edgealpha',0,'facealpha',.3)
contour3(xx,yy,zz0,'linewidth',2);
colormap cool
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% bichromatic elipse
sig2 = sig1 * .65;
params = [A sig1 sig2 orthsig exp bl rot];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz0 = reshape(z,size(xx));
zz1 = reshape(z,size(xx));
idx1 = sign(xx-yy)==1 | sign(xx-yy)==0;
zz1(~idx1) = nan;
zz2 = reshape(z,size(xx));
idx2 = sign(xx-yy)==-1 | sign(xx-yy)==0;
zz2(~idx2) = nan;

name = 'pancolor';
figure(26); clf; hold on; box on; axis square;
set(gcf,'units','pixels','pos',[200 200 500 500])
surf(xx,yy,zz0,'edgealpha',0,'facealpha',.3)
contour(xx,yy,zz1,'linewidth',2)
contour(xx,yy,zz2,'linewidth',2,'linestyle','--')
colormap cool
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Horseshoe %%%
params = [A sig1 sig2 orthsig exp bl rot];
z = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'conicsection_xy');
zz0 = reshape(z,size(xx));
zz1 = reshape(z,size(xx));
idx1 = sign(xx-yy)==1 | sign(xx-yy)==0;
zz1(~idx1) = nan;
zz2 = reshape(z,size(xx));
idx2 = sign(xx-yy)==-1 | sign(xx-yy)==0;
zz2(~idx2) = nan;

name = 'horseshoe';
figure(27); clf; hold on; box on; axis equal tight
set(gcf,'units','pixels','pos',[200 200 500 500])
surf(xx,yy,zz0,'edgealpha',0,'facealpha',.3)
contour(xx,yy,zz1,'linewidth',2)
contour(xx,yy,zz0,'linewidth',2,'linestyle','-')
colormap cool
set(gca,'xtick',[],'ytick',[])
set(gcf,'PaperPositionMode','auto')
export_fig(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-jpg')



%% Figure 3 %%
% Example cell white noise STA, PSTH with rasters, and bubble plot surface

%idx = 78; % sub 1 of simple cell
%idx = 79; % sub 2 of simple cell

%idx = 153; % -lum unichrom
%idx = 82; % +lum unichrom but currently excluded
%idx = 33; % red unichrom
%idx = 70; % green unichrom **

%idx = 92; % bichromatic
%idx = 43; %bichromatic **
idx = 134; % ON-lum

datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
cellsel = ['cellselect(' num2str(idx) ',[])'];
GLMSPopGUI_Tuning(cellsel);
tunfig = get(150,'UserData');
conpanel = get(tunfig.conpanel,'userdata');
surfparams = GLMSPopData{idx+1,strcmp(datatypes,'Surface Parameters')};
cellpan = get(tunfig.cellpanel,'userdata');
surfpan = get(tunfig.surfpanel,'userdata');

%% Figure 3A %%%
% White noise single filter STA
% name = [GLMP.datafile 's' num2str(sub) '(#' num2str(idx) ') WN'];
% figure(30); clf;
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off','name',name,'PaperPositionMode','auto')
% STA = copyobj(cellpan.linfil.axes,30);
% set(STA,'pos',[.1 .1 .8 .8],'title',[])
% axes(STA); axis square
% export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')

% White noise 4 channel filter
chans = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
for n = 1:numel(chans)
   
    name = [GLMP.datafile 's' num2str(sub) '(#' num2str(idx) ') WN ' chans{n}];
    figure(30); clf;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off','name',name,'PaperPositionMode','auto')
    STA = copyobj(cellpan.(chans{n}).axes,30);
    set(STA,'pos',[.1 .1 .8 .8],'title',[])
    axes(STA); axis square
    export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

end


%% Figure 3B %%%
% PSTH
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),bins);
PSTH = (PSTH./sum(GLMP.flashL)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = sum(GLMP.flashL);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot hist
%name = [GLMP.datafile 's' num2str(sub) '(#' num2str(idx) ') PSTH'];
name = 'PSTH';
figure(31); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'renderer','paint',...
    'numbertitle','off','name',name)
set(gca,'tickdir','out')

% right axis
yyaxis right
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
ylabel('Mean Firing Rate (sp/s)')

% Left axis
yyaxis left
x = [];
y = [];
thetarhos = cat(2,GLMP.subunit{sub}.theta,GLMP.subunit{sub}.rho);
[uniq,order] = sortrows(thetarhos);
uniq = flipud(uniq);
order = flipud(order);
for n = 1:numel(order)
    tpts = GLMP.subunit{sub}.normspiketimes{order(n)};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end
h = plot(fliplr(x'),fliplr(y'),'k-');
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Color Direction (rad)')
coldirs = linspace(-180,180,5);
tickspacing = linspace(0,max(smoothPSTH),5);
set(gca,'YTick',tickspacing,'yticklabel',coldirs)

%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Figure 3C %%%
% Bubbble plot surface (this is annoying, but bubbles, contours, and
% surface must all be saved sepeartely for Illustrator to recognize it.

%idx = 70; % green unichrom **
bubblemaker(idx)


%% gradient

x = linspace(-pi,pi,100);
Lc = cos(x);
Lc = (Lc - min(Lc))./(max(Lc)-min(Lc));
Mc = sin(x);
Mc = (Mc - min(Mc))./(max(Mc)-min(Mc));
Sc = repmat(.5,size(Mc));

guns = cat(2,Lc',Mc',Sc');
guns = repmat(guns,1,1,10);
guns = permute(guns,[1 3 2]);

figure(1); clf; hold on;
h = image(guns);
set(gca,'xlim',[min(h.XData) max(h.XData)],'ylim',[min(h.YData) max(h.YData)],...
    'xtick',[],'ytick',[])
name = 'gradient';
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% Gradient second attempt
lummaxcc = .9;
colmaxcc = .09;
thetas = shiftdim(linspace(-pi,pi,361));
rhos = ones(size(thetas));

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = lummaxcc * colmaxcc ./ sqrt((colmaxcc .* cos(thetas-pi/4)).^2 ...
    +(lummaxcc .* sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
Scc = zeros(size(Lcc));

% convert to display colors
DN = GLMSPopData{end,strcmp(datatypes,'DN')};
library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
rawdata = nex2stro([char(library) DN.datafile '.nex']);
bkgndrgb = rawdata.sum.exptParams.bkgndrgb';
bkgndlms =  bkgndrgb * (DN.M);
stimlms = [Lcc Mcc Scc] .* repmat(bkgndlms,size(Lcc,1),1) + repmat(bkgndlms,size(Lcc,1),1);
stimrgb = stimlms * inv(DN.M);

stimrgb = repmat(stimrgb,1,1,10);
stimrgb = permute(stimrgb,[1 3 2]);

figure(2); clf; hold on;
h = image(stimrgb);
set(gca,'xlim',[min(h.XData) max(h.XData)],'ylim',[min(h.YData) max(h.YData)],...
    'xtick',[],'ytick',[])
name = 'gradient';
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Figure 4
% Distribution of normLL vals

%GLMSPopGUI_Params

% Load figure variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

thresh = str2double(conpanel.uicontrols.ndthresh.String);
L = ~poppanel.excludeL;
diffnormLL = poppanel.twoD.bichrom.normLLs(L) - poppanel.oneD.bichrom.normLLs(L);

% populate normLL hist
figure(1); cla; hold on; box on;
maxval = max(diffnormLL);
range = 0:.005:.3;
h = histogram(diffnormLL,range,'edgecolor','w','facecolor','k');
%h = histogram(diffnormLL,thresh:.005:maxval,'edgecolor','w','facecolor','k');
%h = histogram(diffnormLL,0:.005:thresh,'edgecolor','w','facecolor','k');
ylim([0 max(h.Values)*1.1])

%xlim([0 (rndofferr(maxval,2))])
xlim([0 .26])
plot([thresh thresh],[0 max(h.Values)*1.1],'r--')
set(gca,'TickDir','out')

name = 'normLLhist';
set(1,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Figure 5A %%
% 1D polar PD histogram, smoothed polar hist, color modes, color wedges

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

L = poppanel.oneDL;
angs = poppanel.tuning(L);
%angs = poppanel.oneD.bichrom.parvals(L,end-1);
%cardirs = conpanel.fitcard.oneD.means;
angL = conpanel.fitcard.oneD.angL;

disp([num2str(sum(L)) '/' num2str(sum(~poppanel.excludeL)) ' neurons are 1D'])

% Create histogram and smooth
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 45.001*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH = hist(angs,bins);
smoothPSTH = conv(PSTH,gaussfilter,'same');
smoothPSTH = smoothPSTH./max(smoothPSTH) .* max(PSTH); % normalizing by max height of PSTH for plotting purposes
histangs = linspace(-pi,pi,21);

%first derivative
A = spline(bins,smoothPSTH);
M = diag(3:-1:1,1);
D=A;
D.coefs = D.coefs*M;
peaks = nan(4,1);
for n = 1:4
    peaks(n) = fzero(@(bins) ppval(D,bins),conpanel.cardirs(n));
end

figure(31); clf; 
h = polarhistogram(angs,histangs,'Facecolor',[.5 .5 .5],'edgecolor','k'); hold on;
histmax = max(h.Values);
smoothPSTH = smoothPSTH./max(smoothPSTH)*histmax;
h = polarplot(bins,smoothPSTH,'-'); 
set(h,'LineWidth',2,'color',conpanel.oneDcol);

radpm = 22.5/180*pi;
for n = 1:4
    h = polarplot([peaks(n) peaks(n)],[0 histmax]);
    set(h,'color',conpanel.carcols(n,:),'linewidth',3)
    h = polarplot([0 peaks(n)+radpm peaks(n)-radpm 0],...
        [0 histmax*1.1 histmax*1.1 0]);
    set(h,'color',conpanel.carcols(n,:));
end

a = sum(L & any(angL,2));
perc = rndofferr(sum(a)/sum(L)*100,2);
disp([num2str(a) '/' num2str(sum(L)) ' (' num2str(perc)...
    ') neurons are fit by cardinal modes'])

rlim([0 histmax+2])
name = 'bi';
set(31,'PaperPositionMode','auto','renderer','paint')
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Example surfaces %%%

idx = 33; % red 
%idx = 70; % green ** (3C example)
%idx = 43; % -lum
%idx = 134; % +lum

bubblemaker(idx);


%% Figure 5B
% Preferred directions of neurons tested with RS distribution

GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
datatypes = GLMSPopData(1,:);

% parmeter bounds
modmin = 10; % in sp/s
kappamax = 4;

% Modulation
modvals = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    modvals(n-1) = max(GLMP.subunit{sub}.meannspikes)...
        ./ mean(GLMP.subunit{sub}.stimDur) - mean(GLMP.subunit{sub}.blfr);
end

% Kappa
kappavals = nan(size(GLMSPopData,1)-1,1);
kappavals(poppanel.oneD.L) = poppanel.oneD.params(poppanel.oneD.L,end);
kappavals(~poppanel.oneD.L) = poppanel.twoD.params(~poppanel.oneD.L,end);

% datafiles that don't meet inclusion criteria
critL = modvals < modmin | kappavals > kappamax;
disp([num2str(sum(critL)) ' subunits do not meet inclusion criteria.'])

% For selecting one subunit per datafile
filenames = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[~,fnidx,ufnidx] = unique(filenames); 
repsL = ones(size(filenames,1),1); 
repsL(fnidx) = 0; %for 1 sub/cell
repeatsidx = find(repsL);
uniqL = ones(size(filenames,1),1);
for n = 1:numel(repeatsidx)
    idx = find(ufnidx == ufnidx(repeatsidx(n)));
    glmp1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'GLMP')};
    sub1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'Subunit')};
    glmp2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'GLMP')};
    sub2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'Subunit')};
    maxnsp1 = max(glmp1.subunit{sub1}.meannspikes);
    maxnsp2 = max(glmp2.subunit{sub2}.meannspikes);
    if any(critL(idx))
        uniqL(idx) = ~critL(idx);
    elseif maxnsp1 > maxnsp2
        uniqL(idx(2)) = 0;
    else
        uniqL(idx(1)) = 0;
    end
end
disp([num2str(sum(~uniqL & ~critL)) ' secondary subunits excluded from analysis.'])

% Which don't meet criteria?
excludeL = critL | ~uniqL;
disp([num2str(sum(excludeL)) '/' num2str(numel(excludeL)) ' total datafiles excluded from analysis.'])

% RS
RSL = zeros(size(modvals)); 
RSL(1:28) = 1; % hard coded, last RS dataset

% which to include
ndL = poppanel.oneD.L;
angs = poppanel.tuning(ndL & ~RSL & ~excludeL);

% Create histogram and smooth
histangs = linspace(-pi,pi,21);
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 45.001*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH = hist(angs,bins);
smoothPSTH = conv(PSTH,gaussfilter,'same');
smoothPSTH = smoothPSTH./max(smoothPSTH) .* max(PSTH); % normalizing by max height of PSTH for plotting purposes

%first derivative
A = spline(bins,smoothPSTH);
M = diag(3:-1:1,1);
D=A;
D.coefs = D.coefs*M;
peaks = nan(4,1);
for n = 1:4
    peaks(n) = fzero(@(bins) ppval(D,bins),conpanel.cardirs(n));
end

% Plot 1D distribution in background
figure(31); clf; 
polarhistogram(angs,histangs,'displaystyle','stairs','edgecolor','k'); hold on;

% Plot RS distribution in the foreground
angs = poppanel.tuning(ndL & ~excludeL & RSL);
h = polarhistogram(angs,histangs,'facecolor','red','edgecolor','k');
disp([num2str(sum(RSL & ndL & ~excludeL)) '/' num2str(sum(RSL & ~excludeL)) ' RS datasets are 1D'])
histmax = max(h.Values);

% Wedges
radpm = 22.5/180*pi;
for n = 1:4
    h = polarplot([peaks(n) peaks(n)],[0 histmax]);
    set(h,'color',conpanel.carcols(n,:),'linewidth',3)
    h = polarplot([0 peaks(n)+radpm peaks(n)-radpm 0],...
        [0 histmax*1.1 histmax*1.1 0]);
    set(h,'color',conpanel.carcols(n,:));
end

% tidying and saving figure
rlim([0 histmax+2])
name = 'Figure 5C';
set(31,'PaperPositionMode','auto','renderer','paint')
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% 5B Example surfaces %%%

%idx = 8; % green **
%idx = 14; % green
%idx = 24; % red 
idx = 1; % red **

% Unused
%idx = 2; % red 
%idx = 6; % green (sub 2 of DO cell)
%idx = 16; % green
%idx = 17; % bipolar red
%idx = 20; % red
%idx = 22; % bipolar red
%idx = 23; % bipolar red **
%idx = 26; % red


bubblemaker(idx)

%% Figure 5C

% The results of an analysis that measures variance on the fit vs real 
% preferred direction of a 1D model neuron with negative binomial variance.

library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
load([library 'VTM_conicsect_nbinvar'])

% Organize and express in degrees
angs = data.angs./pi*180;
angdiffs = data.angdiff./pi*180;
fitmean = mean(angdiffs,2);
fitvar = var(angdiffs,[],2);
fitstderr = stderr(angdiffs')';
xgridticks = -90:30:90;

% Plot
figure(101); cla; hold on; grid on; box on;
shadedErrorBar(angs,fitmean,fitvar,'r-')
plot(angs,angdiffs,'ko')
plot(angs,fitmean,'r-')
set(gca,'xtick',xgridticks)
xlim([min(angs) max(angs)])
xlabel('Modeled PD (deg)')
ylabel('Error in estimated PD (deg)')

name = 'Figure 5C';
set(gcf,'renderer','paint')
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% Figure 5D %%%
% 1D PD and c50 means and standard deviations
cols = {'+lum' '-lum' 'red chrom' 'green chrom'};

%GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.oneDL;
cardirs = conpanel.fitcard.oneD.means;
angL = conpanel.fitcard.oneD.angL;
name = 'Figure 5D';
radpm = 22.5/180*pi;

nunip = sum(poppanel.params(ndL,3)==0);
perc = nunip/sum(ndL);
disp([num2str(nunip) '/' num2str(sum(ndL)) ' (' num2str(perc) ')'...
    ' of 1D neurons are unipolar.' ])

figure(40); clf; hold on;
xlim([-180 180]); ylim([0 1])
set(gca,'xtick',linspace(-180,180,5))

L = ndL & ~any(angL,2) & ~poppanel.excludeL; % all pts not in wedges
h = plot(poppanel.tuning(L)./pi*180,1./poppanel.params(L,2),'ko'); % PDs not in wedges
set(h,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[1 1 1])
h = plot(poppanel.tuning(L)./pi*180,1./poppanel.params(L,3),'ko'); % APDs not in wedges
set(h,'MarkerEdgeColor',[.5 .5 .5])
for n = 1:4
    L = ndL & angL(:,n) & ~poppanel.excludeL;
    angs = poppanel.tuning(L)./pi*180;
    PDsigs = 1./poppanel.params(L,2);
    apdsigs = 1./poppanel.params(L,3);
    xedge(1) = (cardirs(n)+radpm)/pi*180;
    xedge(2) = (cardirs(n)-radpm)/pi*180;
    apdangs = angs;
%     apdangs = angs+180;
%     if any(apdangs>180)
%         apdangs(apdangs>180) = apdangs(apdangs>180)-360;
%     end

    h = plot([xedge(1) xedge xedge(2)],...
        [0 1 1 0]); % boxes
    set(h,'color',conpanel.carcols(n,:)); hold on;
    h = plot(angs,PDsigs,'o'); % PD sigs within wedge
    set(h,'MarkerEdgeColor','w','MarkerFaceColor',conpanel.carcols(n,:));
    h = plot(apdangs,apdsigs,'o'); % APD sigs within wedge
    set(h,'MarkerEdgeColor',conpanel.carcols(n,:),'MarkerFaceColor','none');
    plot([mean(angs) mean(angs)],[mean(PDsigs)-std(PDsigs) mean(PDsigs)+std(PDsigs)],'k-');
    plot([mean(angs)-std(angs) mean(angs)+std(angs)],[mean(PDsigs) mean(PDsigs)],'k-');
    
    nunip = sum(isinf(apdsigs));
    perc = rndofferr(nunip/numel(apdsigs)*100,2);
    disp([(cols{n}) ' neurons are unipolar: ' num2str(nunip) '/'...
        num2str(numel(apdsigs)) ' (' num2str(perc) '%)'])
    
    temp = mean(PDsigs./apdsigs);
    disp([(cols{n}) ' bipolar ratio ' num2str(temp)])
    
end
set(gcf,'PaperPositionMode','auto','renderer','paint')
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

for n = 1:2
    if n == 1
        L = ndL & any(angL(:,1:2),2);
        str = 'lum';
    else
        L = ndL & any(angL(:,3:4),2);
        str = 'col';
    end
    
    %params = poppanel.oneD.bichrom.parvals(L,:);
    params = poppanel.params(L,:);
    angs = params(:,end)./pi*180;
    pdsigs = params(:,2);
    apdsigs = params(:,3);
    
    
    temp = mean(apdsigs./pdsigs);
    disp([(str) ' bipolar ratio ' num2str(temp)])
    
    temp = apdsigs./pdsigs;
    figure;
    histogram(temp,linspace(0,1,100))
    title(str)
    
end


%% Figure 6 %%
% Example cell white noise STA, PSTH with rasters, and bubble plot surface

%idx = 68; % -lum horseshoe
%idx = 63; % green horsehoe
%idx = 44; % red horseshoe **
%idx = 160; % +lum horseshoe

idx = 135;
%idx = 155; % pan color

datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
cellsel = ['cellselect(' num2str(idx) ',[])'];
GLMSPopGUI_Tuning(cellsel);
tunfig = get(150,'UserData');
conpanel = get(tunfig.conpanel,'userdata');
surfparams = GLMSPopData{idx+1,strcmp(datatypes,'Surface Parameters')};
cellpan = get(tunfig.cellpanel,'userdata');

%% Figure 6A %%%
% % White noise single filter STA
% name = [GLMP.datafile 's' num2str(sub) '(#' num2str(idx) ') WN'];
% figure(30); clf;
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'Renderer','painters','name',name,'PaperPositionMode','auto')
% STA = copyobj(cellpan.linfil.axes,30);
% set(STA,'pos',[.1 .1 .8 .8],'title',[])
% axes(STA); axis square
% export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')

% White noise 4 channel filter
chans = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
for n = 1:numel(chans)
   
    figure(30); clf;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off','name',name,'PaperPositionMode','auto')
    STA = copyobj(cellpan.(chans{n}).axes,30);
    set(STA,'pos',[.1 .1 .8 .8],'title',[])
    axes(STA); axis square
    export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' chans{n}],'-depsc')

end

%% Figure 6B %%%
% PSTH
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),bins);
PSTH = (PSTH./sum(GLMP.flashL)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = sum(GLMP.flashL);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot hist
name = 'PSTH';
figure(31); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'Renderer','painters');
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
x = [];
y = [];
thetarhos = cat(2,GLMP.subunit{sub}.theta,GLMP.subunit{sub}.rho);
[uniq,order] = sortrows(thetarhos);
uniq = flipud(uniq);
order = flipud(order);
for n = 1:numel(order)
    tpts = GLMP.subunit{sub}.normspiketimes{order(n)};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end
[ax,ticax,hLine2] = plotyy(x',y',0,0);
set(ticax(:),'color','k')
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)],'tickdir','out')
%title('PSTH')
%xlabel(ax(1),'Time from Stim Onset (ms)')
%ylabel(ax(1),'Mean Firing Rate (sp/s)')
%ylabel(ax(2),'Color Direction (rad)')
coldirs = linspace(-180,180,5);
set(ax(2),'ylim',[min(coldirs) max(coldirs)],'YTick',coldirs)
%ylabel(ax(2),'Angle in LM plane')
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Figure 6C %%%
% Bubbble plot surface (this is annoying, but bubbles, contours, and
% surface must all be saved sepeartely for Illustrator to recognize it.

idx = 135; % pan color * **
%idx = 155; % pan color
bubblemaker(idx)


%% Figure 7 analysis

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

ndL = poppanel.twoDL;
pcL = ndL & poppanel.twoD.bichrom.eliL;
hsL = ndL & poppanel.twoD.unichrom.eliL;

%allpcangs = poppanel.twoD.bichrom.parvals(ndL,end-1);
%pchsangs = poppanel.params(ndL,end-1);
%pcangs = poppanel.params(pcL,end-1);
hsangs = poppanel.params(hsL,end-1);
hspcangs = poppanel.twoD.bichrom.parvals(hsL,end-1);
histangs = linspace(-pi,pi,21);

figure(60); clf;
h = polarhistogram(hsangs,histangs,'FaceColor',[.5 .5 .5],...
    'edgecolor','k','linewidth',2); hold on;
h = polarhistogram(hspcangs,histangs,'FaceColor','none',...
    'edgecolor','b','linewidth',1); hold on;
legend('HS fit','PC fit')

name = 'HS fit vs PC fit';
set(60,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Analysis for greg

pcnll = poppanel.twoD.bichrom.normLLs(ndL);
hsnll = poppanel.twoD.unichrom.normLLs(ndL);
diffnLL = pcnll-hsnll;
apdc50s = 1./poppanel.twoD.bichrom.parvals(ndL,3);

figure(1); clf; hold on;
plot(diffnLL,apdc50s,'ko')
plot([.07 .07],[0 max(apdc50s)],'r--')
xlim([0 max(diffnLL)])
ylim([0 1])
xlabel('NLL diff')
ylabel('APD c50')

name = '2Dnlldiffs';
set(60,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')
%% Figure 7A %%
% 2D polar PD histogram, smoothed polar hist, color modes, color wedges

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

ndL = poppanel.twoDL;
pcL = ndL & poppanel.twoD.bichrom.eliL & ~poppanel.excludeL;

angmat = nan(sum(pcL),4);
sigmat = nan(sum(pcL),4);

% preferred directions
PDangs = poppanel.tuning(pcL);
PDsigs = 1./poppanel.params(pcL,2);
angmat(:,1) = PDangs;
sigmat(:,1) = PDsigs;

% APD
apdsigs = 1./poppanel.params(pcL,3);
apdangs = PDangs+pi;
if any(apdangs>pi)
    apdangs(apdangs>pi) = apdangs(apdangs>pi)-(2*pi);
end
angmat(:,2) = apdangs;
sigmat(:,2) = apdsigs;
%h = plot(apdangs,apdsigs,'ko');

% Orthogonal
osigs = 1./poppanel.params(pcL,4);
oangs1 = PDangs+(pi/2);
if any(oangs1>pi)
    oangs1(oangs1>pi) = oangs1(oangs1>pi)-(2*pi);
end
angmat(:,3) = oangs1;
sigmat(:,3) = osigs;
%h = plot(oangs1,osigs,'ko');

oangs2 = oangs1+pi;
if any(oangs2 > pi)
    oangs2(oangs2>pi) = oangs2(oangs2>pi)-(2*pi);
end
angmat(:,4) = oangs2;
sigmat(:,4) = osigs';
%h = plot(oangs2,osigs,'ko');

% arrange lowest to highest 
[sortangmat,idx] = sort(angmat,2);
sortsigmat = nan(size(sortangmat));
for n = 1:size(sortsigmat,1)
    sortsigmat(n,:) = sigmat(n,idx(n,:));
end
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

% For mode and wedge
A = spline(bins,smoothPSTH);
M = diag(3:-1:1,1);
D=A;
D.coefs = D.coefs*M;
peak = fzero(@(bins) ppval(D,bins),pi/4);
radpm = 22.5/180*pi;

% Plot histogram  
figure(60); clf;
h = polarhistogram(sortangmat,histangs,'FaceColor','none',...
    'edgecolor','b','linewidth',2); hold on;
histmax = max(h.Values);
h = polarplot([peak peak],[0 histmax]); hold on;
set(h,'color',[.5 .5 .5],'linewidth',3)
h = polarplot([0 peak+radpm peak-radpm 0],...
    [0 histmax*1.1 histmax*1.1 0]);
set(h,'color',[.5 .5 .5]);
% for n = 1:size(sortangmat,1)
%     polarplot([0 sortangmat(n,1)],[0 histmax/2],'k-');
%     polarplot([0 sortangmat(n,2)],[0 histmax/2],'k-');
%     polarplot([0 sortangmat(n,3)],[0 histmax/2],'k-');
%     polarplot([0 sortangmat(n,4)],[0 histmax/2],'k-');
% end
normsmoothpsth = (smoothPSTH./max(smoothPSTH))*histmax+1;
polarplot(bins,normsmoothpsth,'k--','linewidth',2)
rlim([0 histmax+1]);
thetalim([0 90]);

name = 'Figure 7A';
set(60,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

L = any(sortangmat >= peak-radpm & sortangmat <= peak+radpm,2);
perc = rndofferr((sum(L)/numel(L)*100),2);
disp([num2str(perc) '% (' num2str(sum(L)) '/' num2str(numel(L)) ...
    ') principle axes within mode.'])


%% Example surfaces %%%
%idx = 68; % -lum horseshoe
%idx = 63; % green horsehoe
%idx = 44; % red horseshoe **
%idx = 157; % +lum horseshoe
%%idx = 159; % -lum horseshoe

% Pan color examples
%idx = 155; % red pan color (previous
%idx = 52; % -lum pan color *
%idx = 104 % red pc
idx = 135 % red pc * **

bubblemaker(idx);


%% Figure 7B %%%
% 2D PD and c50 means and standard deviations

cols = {'+lum' '-lum' 'red chrom' 'green chrom'};
%GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.twoDL;
cardirs = conpanel.fitcard.oneD.means;
radpm = 22.5/180*pi;
name = 'Figure 7B';

% Pull out pan color cells
pcL = ndL & poppanel.twoD.bichrom.eliL;
hsL = ndL & poppanel.twoD.unichrom.eliL;

% assign variables
prinaxmat = nan(sum(pcL),4);
prinsigmat = nan(sum(pcL),4);

% preferred directions
PDangs = poppanel.tuning(pcL)./pi*180;
PDsigs = 1./poppanel.params(pcL,2);
prinaxmat(:,1) = PDangs;
prinsigmat(:,1) = PDsigs;

% APD
apdsigs = 1./poppanel.params(pcL,3);
apdangs = PDangs+180;
if any(apdangs>180)
    apdangs(apdangs>180) = apdangs(apdangs>180)-360;
end
prinaxmat(:,2) = apdangs;
prinsigmat(:,2) = apdsigs;

% Orthogonal
osigs = 1./poppanel.params(pcL,4);
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
%h = plot(oangs2,osigs,'ko');

[sortangmat,idx] = sort(prinaxmat,2);
sortsigmat = nan(size(sortangmat));
for n = 1:size(sortsigmat,1)
    sortsigmat(n,:) = prinsigmat(n,idx(n,:));
end

% Add padding to both sides (to emulate cyclical plotting)
paddedangmat = cat(2,sortangmat(:,1)-90,sortangmat,sortangmat(:,4)+90);
paddedsigmat = cat(2,sortsigmat(:,4),sortsigmat,sortsigmat(:,1));

% Connect with curved iso-c50 surface
linex = linspace(-pi,pi,129);
sigmat = nan(sum(pcL),numel(linex));
%thetamat = nan(sum(pcL),numel(linex));
pcidx = find(pcL);
for n = 1:numel(pcidx)
    params = poppanel.params(pcidx(n),:);
    pdsig = 1/params(2);
    apdsig = 1/params(3);
    osig = 1/params(4);
    pd = params(end-1);
    shiftedx = linex-pd;
    shiftedx(shiftedx>pi) = shiftedx(shiftedx>pi)-(2*pi);
    shiftedx(shiftedx<-pi) = shiftedx(shiftedx<-pi) + (2*pi);
    posL = shiftedx<pi/2 & shiftedx>-pi/2;
    
    sigmat(n,posL) = (pdsig*osig) ./ sqrt((pdsig*sin(shiftedx(posL))).^2 + (osig*cos(shiftedx(posL))).^2);
    sigmat(n,~posL) = (apdsig*osig) ./ sqrt((apdsig*sin(shiftedx(~posL))).^2 + (osig*cos(shiftedx(~posL))).^2);
        
end


% Figure 
figure(63); clf; hold on; box on;
set(gca,'xtick',linspace(-180,180,7),'tickdir','out')
for n = 1:4
    xedge(1) = (cardirs(n)+radpm)/pi*180;
    xedge(2) = (cardirs(n)-radpm)/pi*180;
    h = plot([xedge(1) xedge xedge(2)],[0 1 1 0]); % boxes
    set(h,'color',conpanel.carcols(n,:)); hold on;
end



shadedErrorBar(linex./pi*180,mean(sigmat),std(sigmat),'b')
alpha(.4)
plot(linex./pi*180,sigmat,'k-')

%idx = 135 % red pc * **
%idxL = find(pcL)==idx;
%plot(linex./pi*180,sigmat(idxL,:),'r-')
%plot(prinaxmat,prinsigmat,'ko')
%h = plot(prinaxmat(~pcL,:),prinsigmat(~pcL,:),'ko');
%set(h,'color',[.5 .5 .5]);
% Plot std along each cardinal axis
% for n = 1:4
%     sigs = sigmat(:,linex==cardaxes(n));
%     h = plot([cardaxes(n)/pi*180 cardaxes(n)/pi*180],[mean(sigs)-std(sigs) mean(sigs)+std(sigs)],'k-');
%     set(h,'linewidth',2)
%     h = plot(cardaxes(n)/pi*180,mean(sigs),'xk');
%     set(h,'linewidth',2)
% %     plot([mean(angs)-std(angs) mean(angs)+std(angs)],[mean(sigs) mean(sigs)],'k-');
% end

xlim([-180 180]); ylim([0 .6])
%xlim([-180 180]); ylim([0 max(prinsigmat(:))*1.1])
title('Pan color cells: axes and c50s')
set(gcf,'PaperPositionMode','auto','numbertitle','off','name',name);
%print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')



%% Figure 8

%idx = 63; % green horsehoe
idx = 44; % red horseshoe **
%idx = 157; % +lum horseshoe
%%idx = 159; % -lum horseshoe

datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
cellsel = ['cellselect(' num2str(idx) ',[])'];
GLMSPopGUI_Tuning(cellsel);
tunfig = get(150,'UserData');
conpanel = get(tunfig.conpanel,'userdata');
surfparams = GLMSPopData{idx+1,strcmp(datatypes,'Surface Parameters')};
cellpan = get(tunfig.cellpanel,'userdata');

%% Figure 8A %%%
% % White noise single filter STA
% name = [GLMP.datafile 's' num2str(sub) '(#' num2str(idx) ') WN'];
% figure(30); clf;
% set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%     'Renderer','painters','name',name,'PaperPositionMode','auto')
% STA = copyobj(cellpan.linfil.axes,30);
% set(STA,'pos',[.1 .1 .8 .8],'title',[])
% axes(STA); axis square
% export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')

% White noise 4 channel filter
chans = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
for n = 1:numel(chans)
   
    name = ['Figure 8 STA ' chans{n}];
    figure(70); clf;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off','name',name,'PaperPositionMode','auto')
    STA = copyobj(cellpan.(chans{n}).axes,70);
    set(STA,'pos',[.1 .1 .8 .8],'title',[])
    axes(STA); axis square
    %export_fig(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')
    print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')
    
end

%% Figure 8B %%%
% PSTH
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.subunit{sub}.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.subunit{sub}.Lcc)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = numel(GLMP.subunit{sub}.Lcc);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot hist
name = 'Figure 8B';
figure(71); clf; hold on; box on;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'Renderer','painters','name',name);
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
x = [];
y = [];
thetarhos = cat(2,GLMP.subunit{sub}.theta,GLMP.subunit{sub}.rho);
[uniq,order] = sortrows(thetarhos);
uniq = flipud(uniq);
order = flipud(order);
for n = 1:numel(order)
    tpts = GLMP.subunit{sub}.normspiketimes{order(n)};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end
[ax,hLine1,hLine2] = plotyy(x',y',0,0);
%hLine1.Color{:} = [1 1 1];
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)],'tickdir','out')
coldirs = linspace(-180,180,13);
set(ax(2),'ylim',[min(coldirs) max(coldirs)],'YTick',coldirs)
ylabel(ax(2),'Angle in LM plane')
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')

%% Figure 8C %%%

idx = 44; % red horseshoe **
bubblemaker(idx)



%% Figure 9A %%
% 2D polar PD histogram, smoothed polar hist, color modes, color wedges

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
name = 'Figure 9A';
ndL = poppanel.twoDL;
hsL = ndL & poppanel.twoD.unichrom.eliL;

% preferred directions
PDangs = poppanel.tuning(hsL);
%PDangs = poppanel.tuning(pcL)+pi;
%PDangs(PDangs>pi) = PDangs(PDangs>pi)-(2*pi);

% Create histogram and smooth
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 60*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH = hist(PDangs,bins);
smoothPSTH = conv(PSTH,gaussfilter,'same');
smoothPSTH = smoothPSTH./max(smoothPSTH) .* max(PSTH); % normalizing by max height of PSTH for plotting purposes
histangs = linspace(-pi,pi,21);

% Calculate color modes
A = spline(bins,smoothPSTH);
M = diag(3:-1:1,1);
D=A;
D.coefs = D.coefs*M;
peaks = nan(4,1);
for n = 1:4
    peaks(n) = fzero(@(bins) ppval(D,bins),conpanel.cardirs(n));
end

% Plot histogram  
figure(81); clf;
h = polarhistogram(PDangs,histangs,'FaceColor','none',...
    'edgecolor','k','linewidth',2); hold on;
histmax = max(h.Values);

% plot color modes and wedges
radpm = 22.5/180*pi;
for n = 1:4
    h = polarplot([peaks(n) peaks(n)],[0 histmax]); hold on;
    set(h,'color',conpanel.carcols(n,:),'linewidth',3)
    h = polarplot([0 peaks(n)+radpm peaks(n)-radpm 0],...
        [0 histmax*1.1 histmax*1.1 0]);
    set(h,'color',conpanel.carcols(n,:));
end

% grey wedge, from PA fit (8a)
% radpm = 22.5/180*pi;
% for n = 1:4
%     p = peak+(pi/2*n);
%     h = polarplot([0 p+radpm p-radpm 0],...
%         [0 histmax*1.1 histmax*1.1 0]);
%     set(h,'color',[.5 .5 .5]);
% end

normsmoothpsth = (smoothPSTH./max(smoothPSTH))*histmax*1.1;
polarplot([bins bins(1)],[normsmoothpsth normsmoothpsth(1)],'k--','linewidth',1)
rlim([0 histmax+1]);
%thetalim([0 90]);

set(gcf,'PaperPositionMode','auto','renderer','paint',...
    'numbertitle','off','name',name)
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% L = any(angL,2);
% perc = rndofferr((sum(L)/numel(L)*100),2);
% disp([num2str(perc) '% (' num2str(sum(L)) '/' num2str(numel(L)) ...
%    ') principle axes within mode.'])


%% Example surfaces %%%

%idx = 63; % green horsehoe
%idx = 44; % red horseshoe **
%idx = 158; % +lum horseshsoe
idx = 68; % -lum horseshoe

bubblemaker(idx);


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
% angL = nan(sum(hsL),4);
% for n = 1:4
%     xedge(1) = (peaks(n)+radpm)/pi*180;
%     xedge(2) = (peaks(n)-radpm)/pi*180;
%     angL(:,n) = PDangs < xedge(1) & PDangs > xedge(2);
%     h = plot([xedge(1) xedge xedge(2) xedge(1)],[0 1 1 0 0]); % boxes
%     set(h,'color',conpanel.carcols(n,:)); hold on;
% end

% Calculate iso-c50 contours
%linex = linspace(-pi,pi,181);
linex = linspace(-pi,pi,9)
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
% x1 = nanmean(sigmat);
% x2 = nanstd(sigmat);
% order = 3;
% framelen = 31;
% smoothmean = sgolayfilt(x1,order,framelen);
% smoothstd = sgolayfilt(x2,order,framelen);

% Plot smoothed
h = plot(linex./pi*180,sigmat,'k-');
set(h,'color',[.85 .85 .85])
shadedErrorBar(linex./pi*180,smoothmean,smoothstd,'b')
%plot(linex./pi*180,nanmean(sigmat),'r')
alpha(.4)

% wedge points
L = any(angL,2);
plot(PDangs(L),PDsigs(L),'ko','markerfacecolor','k')
oangs = PDangs+90;
oangs(oangs>180) = oangs(oangs>180)-360;
plot(oangs(L),osigs(L),'ko')
oangs = PDangs-90;
oangs(oangs<-180) = oangs(oangs<-180)+360;
plot(oangs(L),osigs(L),'ko')

% ungrouped pts
L = ~any(angL,2);
plot(PDangs(L),PDsigs(L),'o','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
oangs = PDangs+90;
oangs(oangs>180) = oangs(oangs>180)-360;
plot(oangs(L),osigs(L),'o','color',[.5 .5 .5])
oangs = PDangs-90;
oangs(oangs<-180) = oangs(oangs<-180)+360;
plot(oangs(L),osigs(L),'o','color',[.5 .5 .5])

% Figure housekeeping
xlim([-180 180]); ylim([0 1])
title('Horseshoe cells: principle axes and c50s')
set(gcf,'PaperPositionMode','auto','numbertitle','off','name',name);
print(['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%% For figure, use 9 point linex
lummean = nanmean(sigmat(:,[2,6]),2)
colmean = nanmean(sigmat(:,[4,8]),2)
mean(lummean)
mean(colmean)
pval = signrank(colmean,lummean)

%% Discussion Figure?
% Several 1D neurons making a horseshoe.
% Horseshoe + 1d = pancolor.

%idx = 70; % green chrom
%idx = 77; % -lum
%idx = 111; % +lum
%idx = 63; % green horseshoe
%idx = 130; % red chrom
idx = 155; % pan color
idx = 52; %

GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
cellsel = ['cellselect(' num2str(idx) ',[])'];
%GLMSPopGUI_Tuning(cellsel);
tunfig = get(150,'UserData');
conpanel = get(tunfig.conpanel,'userdata');
surfparams = GLMSPopData{idx+1,strcmp(datatypes,'Surface Parameters')};
cellpan = get(tunfig.cellpanel,'userdata');

figure(33); clf; % hold on; box on;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
colormap('cool')
params = conpanel.oneD.params(conpanel.selectedidx,:);
meannsp = GLMP.subunit{sub}.meannspikes;
ticks = -.6:.2:.6;

% Define thetas and rhos for the whole grid
thetas = linspace(-pi,pi,361)';
majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.05;
minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.5;
nom = majorax * minorax;
denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
rhos = nom ./ sqrt(denom);
scalars = 0:.05:1;
rhosgrid = rhos * scalars;
thetasgrid = repmat(thetas,1,numel(scalars));
[x,y] = pol2cart(thetasgrid,rhosgrid);
surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surfparams.twoD.surftype);
surface = reshape(surface,size(rhosgrid));
axlim = max(x(:));

% (1) Surface
% h = surf(x,y,surface); hold on;
% set(h,'edgecolor','none');
% alpha(.3);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 

% (2) Bubbles 
polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
    mn = meannsp(i)/max(meannsp)*50+5; hold on;
    h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
    set(h,'MarkerFaceColor','none','MarkerSize',mn)
end
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 

% % (3) contours and pd
% contrhos = linspace(0,max(GLMP.subunit{sub}.rho),50);
% [tt,rr] = meshgrid(thetas,contrhos);
% [xx,yy] = pol2cart(tt,rr);
% contsurf = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],surfparams.twoD.surftype);
% contsurf = reshape(contsurf,size(rr));
% rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
% PD = rotmat * [0 .4]';
% contour3(xx,yy,contsurf,5,'linewidth',2);
% contour3(x,y,surface,5,'linewidth',2); hold on;
% plot([0 PD(1)],[0 PD(2)],'k-')
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

% axis stuff
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
name = 'DiscussionSurf';
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/Example Cells/' name],'-depsc')




%% Figure 10 RF Structures

% 1D match
%idx = 111; % primary 
idx = 112; % secondary

% 2D match
%idx = 149; % primary
%idx = 150; % secondary

% Mismatch
idx = 80; % primary
%idx = 81; % secondary

DN = GLMSPopData{idx+1,strcmp(datatypes,'DN')};
%filters = GLMSPopData{idx+1,strcmp(datatypes,'DN Filters')};
frames = 6:9; % which frames to average over
nstix = DN.NStixGrid(1);

% +L+M
stas = DN.stats.pLpM.STA(:,frames);
stas = reshape(stas,[nstix nstix 3 numel(frames)]);
stas = stas(:,:,1,:);
sta = mean(stas,4);
pLpMsta = repmat(sta,[1 1 3]);

% -L-M
stas = DN.stats.mLmM.STA(:,frames);
stas = reshape(stas,[nstix nstix 3 numel(frames)]);
stas = stas(:,:,1,:);
sta = mean(stas,4);
mLmMsta = abs(repmat(sta,[1 1 3]));

% +L-M
stas = DN.stats.pLmM.STA(:,frames);
stas = reshape(stas,[nstix nstix 3 numel(frames)]);
stas = stas(:,:,1,:);
sta = mean(stas,4);
pLmMsta = abs(repmat(sta,[1 1 3]));

% -L+M
stas = DN.stats.mLpM.STA(:,frames);
stas = reshape(stas,[nstix nstix 3 numel(frames)]);
stas = stas(:,:,1,:);
sta = mean(stas,4);
mLpMsta = abs(repmat(sta,[1 1 3]));


minval = min([min(pLpMsta(:)) min(mLmMsta(:)) min(pLmMsta(:)) min(mLpMsta(:))]);
maxval = max([max(pLpMsta(:)) max(mLmMsta(:)) max(pLmMsta(:)) max(mLpMsta(:))]);


figure(1); clf;
normsta = (pLpMsta-minval) ./ (maxval-minval);
image(normsta);
axis('square');
set(gca,'xtick',[],'ytick',[])
set(1,'PaperPositionMode','auto')
%export_fig('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/pLpM','-eps')

figure(2); clf;
normsta = (mLmMsta-minval) ./ (maxval-minval);
image(normsta);
axis('square');
set(gca,'xtick',[],'ytick',[])
set(1,'PaperPositionMode','auto')
%export_fig('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/mLmM','-eps')

figure(3); clf;
normsta = (pLmMsta-minval) ./ (maxval-minval);
image(normsta);
axis('square');
set(gca,'xtick',[],'ytick',[])
set(1,'PaperPositionMode','auto')
export_fig('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/pLmM','-eps')

figure(4); clf;
normsta = (mLpMsta-minval) ./ (maxval-minval);
image(normsta);
axis('square');
set(gca,'xtick',[],'ytick',[])
set(1,'PaperPositionMode','auto')
export_fig('/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/mLpM','-eps')

%% Figure 10 surfaces
% sub1 and sub2 agree/disagree


% 1D match
%idx = 111; % primary 
%idx = 112; % secondary

% 2D match
%idx = 149; % primary
%idx = 150; % secondary

% Mismatch
%idx = 80; % primary
idx = 81; % secondary

bubblemaker(idx)

%%

name = 'Figure 10 Surf';

GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};

ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
params = poppanel.params(idx,:);

% Bubbble plot surface (this is annoying, but bubbles, contours, and
% surface must all be saved sepeartely for Illustrator to recognize it.
meannsp = GLMP.subunit{sub}.meannspikes;
ticks = -.8:.2:.8;
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
figure(33); clf; % hold on; box on;
colormap(cmap)
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
h = surf(x,y,surface); hold on;
set(h,'edgecolor','none');
alpha(.3);
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

name = 'surf';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%%% (2) contours and pd %%%
figure(2); clf; colormap(cmap)
rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
PD = rotmat * [0 .5]';
contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
plot([0 PD(1)],[0 PD(2)],'k-')
set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

% save fig
name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

% (3) Bubbles
figure(3); clf;
polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
scalar = 30;
maxnsp = max(cat(1,GLMP.subunit{1}.meannspikes,GLMP.subunit{2}.meannspikes));
for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
    mn = meannsp(i)/maxnsp*scalar+(scalar/10);
    h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko'); 
    set(h,'MarkerFaceColor','none','MarkerSize',mn)
end
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
L = GLMP.subunit{sub}.meannspikes == max(GLMP.subunit{sub}.meannspikes);
h = polar(GLMP.subunit{sub}.uniquetheta(L),GLMP.subunit{sub}.uniquerho(L),'ko');
set(h,'MarkerFaceColor','none','MarkerSize',meannsp(L)/maxnsp * scalar + (scalar/10))
L = GLMP.subunit{sub}.meannspikes == min(GLMP.subunit{sub}.meannspikes);
k = polar(GLMP.subunit{sub}.uniquetheta(L),GLMP.subunit{sub}.uniquerho(L),'ko');
set(k,'MarkerFaceColor','none','MarkerSize',scalar/10)
frs = num2str(cat(1,round(max(GLMP.subunit{sub}.meanfr)),round(min(GLMP.subunit{sub}.meanfr))));
legend([h;k],frs,'location','northwest')

% save fig
name = 'bubbles';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
%print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

%% Figure 10 (unused)
% PD histogram of secondary subs

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

excludeL = ~poppanel.modL | ~poppanel.kappaL | poppanel.RSL;
filenames = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
subs = [GLMSPopData{2:end,strcmp(datatypes,'Subunit')}]';
secsubL = zeros(size(filenames,1),1);

% only look at datafiles with 2 subs that meet criteria
nonexfilenames = filenames(~excludeL);
nonexsubs = subs(~excludeL);

% If both subs meet criteria, take lesser responsive one
[~,ufnidx,fnidx] = unique(nonexfilenames); 
for n = 1:numel(ufnidx)
    idx = find(ufnidx(n) == fnidx);
    if numel(idx) > 1
        popidx = find(strcmp(nonexfilenames{idx(1)},filenames))+1;
        glmp1 = GLMSPopData{popidx(1),strcmp(datatypes,'GLMP')};
        sub1 = GLMSPopData{popidx(1),strcmp(datatypes,'Subunit')};
        glmp2 = GLMSPopData{popidx(2),strcmp(datatypes,'GLMP')};
        sub2 = GLMSPopData{popidx(2),strcmp(datatypes,'Subunit')};
        maxnsp1 = max(glmp1.subunit{sub1}.meannspikes);
        maxnsp2 = max(glmp2.subunit{sub2}.meannspikes);
        if maxnsp1 > maxnsp2
            secsubL(popidx(2)-1) = 1;
        else
            secsubL(popidx(1)-1) = 1;
        end
    end
end
disp([num2str(sum(secsubL)) ' secondary subunits meet criteria.'])


% 1D Preferred Directions
ndL = poppanel.oneD.L;
L = secsubL & ndL;
disp([num2str(sum(L)) '/' num2str(sum(secsubL)) ' neurons are 1D'])
angs = poppanel.tuning(L);
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 45.001*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH = hist(angs,bins);
smoothPSTH = conv(PSTH,gaussfilter,'same');
smoothPSTH = smoothPSTH./max(smoothPSTH) .* max(PSTH); % normalizing by max height of PSTH for plotting purposes

% Set up figure
figure(31); clf;
set(gcf,'numbertitle','off','Name','Secondary Subunit PDs')
histangs = linspace(-pi,pi,21);

% Plot 1D distribution in background
angs = poppanel.tuning(ndL & ~poppanel.excludeL);
h = polarhistogram(angs,histangs,'displaystyle','stairs','edgecolor','k'); hold on;
histmax = max(h.Values);
set(gca,'RLim',[0 histmax+1]);

% Wedges
radpm = 22.5/180*pi;
peaks = conpanel.fitcard.oneD.means;
for n = 1:4
    h = polarplot([peaks(n) peaks(n)],[0 histmax]);
    set(h,'color',conpanel.carcols(n,:),'linewidth',3)
    h = polarplot([0 peaks(n)+radpm peaks(n)-radpm 0],...
        [0 histmax+1 histmax+1 0]);
    set(h,'color',conpanel.carcols(n,:));
end

% Plot secondary subunits
angs = poppanel.tuning(L);
h = polarhistogram(angs,histangs,'facecolor',conpanel.oneDcol,'edgecolor','k');
histmax = max(h.Values);
% smoothPSTH = smoothPSTH./max(smoothPSTH)*histmax;
% h = polarplot(bins,smoothPSTH,'-'); 
% set(h,'LineWidth',2,'color',conpanel.oneDcol);


% 2D Principal Axes
L = secsubL & ~poppanel.oneD.L;
disp([num2str(sum(L)) '/' num2str(sum(secsubL)) ' neurons are 2D'])


fns = cat(1,GLMSPopData{find(secsubL)+1,1});

sum(fns(:,1)=='N')
sum(fns(:,1)=='M')




%% Figure 11 A
% Distribution of flash/Gabor correlation with two examples, high and low

GLMSPopGUI_Gabors
gabfig = get(16,'userdata');
conpanel = get(gabfig.conpanel,'userdata');

figure(3); clf; hold on; box on;
histogram(conpanel.corranal.corrcoef,-1:.05:1);
set(gca,'xlim',[.5 1],'Tickdir','out')

set(gcf,'PaperPositionMode','auto','renderer','paint')

name = 'Figure 11 hist';
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

maxresptab = nan(numel(conpanel.gaboridx),2);
for n = 1:numel(conpanel.gaboridx)
    maxresptab(n,1) = max(conpanel.corranal.flash.allresps{n});
    maxresptab(n,2) = max(conpanel.corranal.gabor.allresps{n});
end
L = conpanel.corranal.corrcoef > .8;
min(maxresptab(L,:))

%% Figure 11 B
% Bubble plots

GLMSPopGUI_Gabors
gabfig = get(16,'userdata');
conpanel = get(gabfig.conpanel,'userdata');

%idx = 84; % Highest correlation
idx = 148; % High correlation

%idx = 97; % low correlation


gabidx = find(conpanel.gaboridx==idx);
stim = conpanel.corranal.stim{gabidx};
flashresps = conpanel.corranal.flash.allresps{gabidx};
gaborresps = conpanel.corranal.gabor.allresps{gabidx};
allresps = cat(1,flashresps,gaborresps);
maxresp = max(allresps);
minresp = min(allresps);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};

% Set up figure
figure(4); clf; 
axlim = max(GLMP.subunit{1}.Lcc);
ticks = -.8:.2:.8;

% (1) Overlaid bubble plot
polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
L = find(allresps == maxresp);
mn = 30+3;
allstim = cat(1,stim,stim);
k = polar(allstim(L,1),allstim(L,2),'ko');
set(k,'MarkerFaceColor','y','MarkerSize',mn)

% Gabors
for i = 1:size(stim,1)
    mn = gaborresps(i)/maxresp*30+3;
    h = polar(stim(i,1),stim(i,2),'o');
    set(h,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','MarkerSize',mn)
end

% Flashes
for i = 1:size(stim,1)
    mn = flashresps(i)/maxresp*30+3;
    h = polar(stim(i,1),stim(i,2),'ko');
    set(h,'MarkerFaceColor','none','MarkerSize',mn)
end

% Axis stuff
axis equal square tight; box on;
set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])
legend(k,['max resp = ' num2str(max(gaborresps))],'location','northwest')

% (2) Axis ticks for box (no surface or contour)
% box on; grid off;
% axis equal square tight
% set(gca,'tickdir','out',...
%     'xtick',ticks,'ytick',ticks,...
%     'xlim',[-axlim axlim],'ylim',[-axlim axlim]);


print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/gabors'],'-depsc')


%% Figure 12: Misclassification

% 1D but looks 2D
idx = 87; % looks pan color trech fit **
%idx = 94; % looks pan color trench fit
%idx = 128; % looks pan color trench fit
%idx = 131; % looks  pan color, trench fit
%idx = 132; % looks horseshoe
%idx = 145; % looks pan color, trench fit *

% 2D but looks 1D
%idx = 32; % looks lum but HS fit **
%idx = 34; % looks lum but HS fit 
%idx = 107; % looks lum but HS fit *

% cone weights
%idx = 42; % cone opponent lum cel **
%idx = 59; % cone opponent lum cell
%idx = 114; % cone opponent lum cell


bubblemaker(idx)


%% Discussion figure: 1D/2D misclassificiations

idx = 87; % nice
idx = 91; % nice
idx = 94; % awesome
idx = 107; % nice
idx = 128; % nice
idx = 131; % looks like pan color, nice
idx = 132; % looks like horseshoe
idx = 138; % weak
idx = 145; % nice




