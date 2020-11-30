function GLMSPopGUI_Params(~,~)

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    library = 'C:\Users\Patty\Dropbox\Patrick\GLMS Data\';
end
%load([library 'GLMSPopData.mat'])

% Main Script
setUpFig();
SetUpHists();
resetplots();

end


%%% Workhorse functions %%%

function PlotTuningHists()

% Load figure variables
ParamsFig = get(gcf,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
%ndparpanel = get(ParamsFig.ndparpanel,'userdata');

% Which monkey
if poppanel.uicontrols.whichMonk.Value == 1
    monkL = ones(size(poppanel.monk.NutL));
elseif poppanel.uicontrols.whichMonk.Value == 2
    monkL = poppanel.monk.NutL;
elseif poppanel.uicontrols.whichMonk.Value == 3
    monkL = poppanel.monk.MauiL;
end

% For plotting 1 subunit per cell
%datafiles = conpanel.table.Data(:,1);
datafiles = poppanel.filenames;
%[uniquefile,a,b] = unique(datafiles); uniqL = zeros(size(datafiles)); uniqL(a) = 1; %for 1 sub/cell
uniqL = ones(size(datafiles)); % all subunits

% Which 1D surface
if poppanel.uicontrols.which1Dsurf.Value == 1
    surf1dL = poppanel.oneD.L;
elseif poppanel.uicontrols.which1Dsurf.Value == 2
    surf1dL = poppanel.oneD.unichrom.L;
elseif poppanel.uicontrols.which1Dsurf.Value == 3
    surf1dL = poppanel.oneD.bichrom.L;
end

% Which 2D Surface
if poppanel.uicontrols.which2Dsurf.Value == 1
    surf2dL = ones(size(poppanel.oneD.L));
elseif poppanel.uicontrols.which2Dsurf.Value == 2
    surf2dL = poppanel.twoD.unichrom.hypL;
elseif poppanel.uicontrols.which2Dsurf.Value == 3
    surf2dL = poppanel.twoD.bichrom.hypL;
elseif poppanel.uicontrols.which2Dsurf.Value == 4
    surf2dL = poppanel.twoD.unichrom.eliL;
elseif poppanel.uicontrols.which2Dsurf.Value == 5
    surf2dL = poppanel.twoD.bichrom.eliL;
end

% Which 1d and 2d surfs to show
poppanel.oneDL = poppanel.oneD.L & monkL & surf1dL & uniqL & ~poppanel.excludeL;
poppanel.twoDL = ~poppanel.oneD.L & monkL & surf2dL & uniqL & ~poppanel.excludeL;

% Bin centers, edges, and counts
bins = linspace(-pi,pi,poppanel.nbins*2+1);
poppanel.hist.binedges = bins(1:2:end);
poppanel.hist.bincenters = bins(2:2:end);
poppanel.hist.counts.oneD = histcounts(poppanel.tuning(poppanel.oneDL),poppanel.hist.binedges);
poppanel.hist.counts.twoD = histcounts(poppanel.tuning(poppanel.twoDL),poppanel.hist.binedges);
radpm = str2double(conpanel.uicontrols.wedgesize.String)/2/180*pi;

% Plot preferred axes of 1D cells
angs = poppanel.tuning(poppanel.oneDL);
axes(poppanel.axes.oneD); cla;
for n = 1:4
    h = polar([0 conpanel.fitcard.oneD.means(n)+radpm conpanel.fitcard.oneD.means(n)-radpm 0],...
        [0 max(poppanel.hist.counts.oneD)*1.1 max(poppanel.hist.counts.oneD)*1.1 0]);
    set(h,'color',conpanel.carcols(n,:)); hold on;
end
poppanel.hist.oneD = rose(angs,poppanel.nbins); hold on;
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13)); % Deleting angle labels and outermost rho label
end
axis tight
set(poppanel.hist.oneD,'color',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs)
title('1D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
% angs = conpanel.tuning(conpanel.oneD.L & monkL & surf1dL & uniqL & ~conpanel.exclude & conpanel.oneD.selectedangs);
% if sum(conpanel.oneD.selectedangs > 0)
%     poppanel.hist.oneDselected = rose2(angs,poppanel.nangs,conpanel.oneDcol,...
%         'edgecolor',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs);
% end

% Plot axes of 2D cells
angs = poppanel.tuning(poppanel.twoDL);
axes(poppanel.axes.twoD); cla;
for n = 1:4
    h = polar([0 conpanel.fitcard.twoD.means(n)+radpm conpanel.fitcard.twoD.means(n)-radpm 0],...
        [0 max(poppanel.hist.counts.twoD)*1.1 max(poppanel.hist.counts.twoD)*1.1 0]);
    set(h,'color',conpanel.carcols(n,:)); hold on;
end
poppanel.hist.twoD = rose(angs,poppanel.nbins); hold on; axis tight
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13));% Deleting angle labels and outermost rho label
end
set(poppanel.hist.twoD,'color',conpanel.twoDcol,'ButtonDownFcn',@SelectTuningDirs)
title('2D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
% angs = conpanel.tuning(~conpanel.oneD.L & monkL & surf2dL & uniqL & ~conpanel.exclude & conpanel.twoD.selectedangs);
% if sum(conpanel.twoD.selectedangs > 0)
%     poppanel.hist.twoDselected = rose2(angs,poppanel.nangs,conpanel.twoDcol);
%     set(poppanel.hist.twoDselected,'ButtonDownFcn',@SelectTuningDirs);
% end


% Save Data
set(ParamsFig.poppanel,'UserData',poppanel)
set(ParamsFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',ParamsFig)

end

function PlotParamHists()

% Load figure variables
ParamsFig = get(gcf,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndparpanel = get(ParamsFig.ndparpanel,'userdata');
onedparpanel = get(ParamsFig.onedparpanel,'userdata');
twodparpanel = get(ParamsFig.twodparpanel,'userdata');

%%% 1D vs 2D %%%

% Upper Asymptote
oneDvals = poppanel.oneD.params(poppanel.oneDL,1);
twoDvals = poppanel.twoD.params(poppanel.twoDL,1);
maxval = max(cat(1,oneDvals,twoDvals));
bins = linspace(0,maxval,2*poppanel.nbins+1);
ndparpanel.upperA.binedges = bins(1:2:end);
ndparpanel.upperA.bincenters = bins(2:2:end);
ndparpanel.upperA.counts.oneD = histcounts(oneDvals,ndparpanel.upperA.binedges);
ndparpanel.upperA.counts.twoD = histcounts(twoDvals,ndparpanel.upperA.binedges);
axes(ndparpanel.upperA.axes); cla; hold on; grid on; box on;
ndparpanel.upperA.hist.oneD = bar(ndparpanel.upperA.bincenters,ndparpanel.upperA.counts.oneD);
set(ndparpanel.upperA.hist.oneD,'FaceColor',conpanel.oneDcol,'EdgeColor',conpanel.oneDcol,...
    'facealpha',0.6)
ndparpanel.upperA.hist.twoD = bar(ndparpanel.upperA.bincenters,ndparpanel.upperA.counts.twoD);
set(ndparpanel.upperA.hist.twoD,'FaceColor',conpanel.twoDcol,'EdgeColor',conpanel.twoDcol,...
    'facealpha',0.6)
xlim([0 max(ndparpanel.upperA.binedges)]);
ylim([0 max(cat(2,ndparpanel.upperA.counts.oneD,ndparpanel.upperA.counts.twoD))]);

% Basline
oneDvals = poppanel.oneD.params(poppanel.oneDL,end-2);% hard coded mapping
twoDvals = poppanel.twoD.params(poppanel.twoDL,end-2);
maxval = max(cat(1,oneDvals,twoDvals));
bins = linspace(0,maxval,2*poppanel.nbins+1);
ndparpanel.bl.binedges = bins(1:2:end);
ndparpanel.bl.bincenters = bins(2:2:end);
ndparpanel.bl.counts.oneD = histcounts(oneDvals,ndparpanel.bl.binedges);
ndparpanel.bl.counts.twoD = histcounts(twoDvals,ndparpanel.bl.binedges);
axes(ndparpanel.bl.axes); cla; hold on; grid on; box on;
h1 = bar(ndparpanel.bl.bincenters,ndparpanel.bl.counts.oneD);
set(h1,'FaceColor',conpanel.oneDcol,'EdgeColor',conpanel.oneDcol,...
    'facealpha',0.6)
h2 = bar(ndparpanel.bl.bincenters,ndparpanel.bl.counts.twoD);
set(h2,'FaceColor',conpanel.twoDcol,'EdgeColor',conpanel.twoDcol,...
    'facealpha',0.6)
xlim([0 max(ndparpanel.bl.binedges)]);
ylim([0 max(cat(2,ndparpanel.bl.counts.oneD,ndparpanel.bl.counts.twoD))]);

% Sigmas (just taking principle axis sigma for every datafile)
oneDvals = 1./poppanel.oneD.params(poppanel.oneDL,2);
twoDvals = 1./poppanel.twoD.params(poppanel.twoDL,2);
maxval = max(cat(1,oneDvals,twoDvals));
bins = linspace(0,maxval,2*poppanel.nbins+1);
ndparpanel.sig.binedges = bins(1:2:end);
ndparpanel.sig.bincenters = bins(2:2:end);
ndparpanel.sig.counts.oneD = histcounts(oneDvals,ndparpanel.sig.binedges);
ndparpanel.sig.counts.twoD = histcounts(twoDvals,ndparpanel.sig.binedges);
axes(ndparpanel.sig.axes); cla; hold on; grid on; box on;
h1 = bar(ndparpanel.sig.bincenters,ndparpanel.sig.counts.oneD);
set(h1,'FaceColor',conpanel.oneDcol,'EdgeColor',conpanel.oneDcol,...
    'facealpha',0.6)
h2 = bar(ndparpanel.sig.bincenters,ndparpanel.sig.counts.twoD);
set(h2,'FaceColor',conpanel.twoDcol,'EdgeColor',conpanel.twoDcol,...
    'facealpha',0.6)
xlim([0 max(ndparpanel.sig.binedges)]);
ylim([0 max(cat(2,ndparpanel.sig.counts.oneD,ndparpanel.sig.counts.twoD))]);

% Exponent
oneDvals = poppanel.oneD.params(poppanel.oneDL,5);
twoDvals = poppanel.twoD.params(poppanel.twoDL,5);
maxval = max(cat(1,oneDvals,twoDvals));
bins = linspace(0,maxval,2*poppanel.nbins+1);
ndparpanel.exp.binedges = bins(1:2:end);
ndparpanel.exp.bincenters = bins(2:2:end);
ndparpanel.exp.counts.oneD = histcounts(oneDvals,ndparpanel.exp.binedges);
ndparpanel.exp.counts.twoD = histcounts(twoDvals,ndparpanel.exp.binedges);
axes(ndparpanel.exp.axes); cla; hold on; grid on; box on;
h1 = bar(ndparpanel.exp.bincenters,ndparpanel.exp.counts.oneD);
set(h1,'FaceColor',conpanel.oneDcol,'EdgeColor',conpanel.oneDcol,...
    'facealpha',0.6)
h2 = bar(ndparpanel.exp.bincenters,ndparpanel.exp.counts.twoD);
set(h2,'FaceColor',conpanel.twoDcol,'EdgeColor',conpanel.twoDcol,...
    'facealpha',0.6)
xlim([0 max(ndparpanel.exp.binedges)]);
ylim([0 max(cat(2,ndparpanel.exp.counts.oneD,ndparpanel.exp.counts.twoD))]);

% Kappa
oneDvals = poppanel.oneD.params(poppanel.oneDL,end);
twoDvals = poppanel.twoD.params(poppanel.twoDL,end);
maxval = max(cat(1,oneDvals,twoDvals));
bins = linspace(0,maxval,2*poppanel.nbins+1);
ndparpanel.kappa.binedges = bins(1:2:end);
ndparpanel.kappa.bincenters = bins(2:2:end);
ndparpanel.kappa.counts.oneD = histcounts(oneDvals,ndparpanel.kappa.binedges);
ndparpanel.kappa.counts.twoD = histcounts(twoDvals,ndparpanel.kappa.binedges);
axes(ndparpanel.kappa.axes); cla; hold on; grid on; box on;
h1 = bar(ndparpanel.kappa.bincenters,ndparpanel.kappa.counts.oneD);
set(h1,'FaceColor',conpanel.oneDcol,'EdgeColor',conpanel.oneDcol,...
    'facealpha',0.6)
h2 = bar(ndparpanel.kappa.bincenters,ndparpanel.kappa.counts.twoD);
set(h2,'FaceColor',conpanel.twoDcol,'EdgeColor',conpanel.twoDcol,...
    'facealpha',0.6)
xlim([0 max(ndparpanel.kappa.binedges)]);
ylim([0 max(cat(2,ndparpanel.kappa.counts.oneD,ndparpanel.kappa.counts.twoD))]);


%%% Color Directions %%%
nds = {'oneD' 'twoD'};

% Plot 1D and 2D params by color direction
for i = 1:2
    whichnd = nds{i};
    if i == 1
        selectedL = poppanel.oneDL;
        cardirs = conpanel.fitcard.oneD.means;
        angL = conpanel.fitcard.oneD.angL;
    else
        selectedL = poppanel.twoDL;
        cardirs = conpanel.fitcard.twoD.means;
        angL = conpanel.fitcard.twoD.angL;
    end
    
    % Upper Asymptote
    if i == 1
        ax = onedparpanel.upperA.axes;
    else
        ax = twodparpanel.upperA.axes;
    end
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    counts = nan(4,poppanel.nbins);
    for n = 1:numel(cardirs)
        L = selectedL & angL(:,n);
        counts(n,:) = histcounts(poppanel.(whichnd).params(L,1),ndparpanel.upperA.binedges);
    end
    h = bar(ax,ndparpanel.upperA.bincenters,counts','stacked');
    for n = 1:4
        set(h(n),'FaceColor',conpanel.carcols(n,:))
    end
    xlim(ax,[0 max(ndparpanel.upperA.binedges)]);
    ylim(ax,[0 max(cat(2,ndparpanel.upperA.counts.oneD,ndparpanel.upperA.counts.twoD))]);
    
    
    % Baseline
    if i == 1
        ax = onedparpanel.bl.axes;
    else
        ax = twodparpanel.bl.axes;
    end
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    counts = nan(4,poppanel.nbins);
    for n = 1:numel(cardirs)
        L = selectedL & angL(:,n);
        counts(n,:) = histcounts(poppanel.(whichnd).params(L,6),ndparpanel.bl.binedges);
    end
    h = bar(ax,ndparpanel.bl.bincenters,counts','stacked');
    for n = 1:4
        set(h(n),'FaceColor',conpanel.carcols(n,:))
    end
    xlim(ax,[0 max(ndparpanel.bl.binedges)]);
    ylim(ax,[0 max(cat(2,ndparpanel.bl.counts.oneD,ndparpanel.bl.counts.twoD))]);
    
    % Sigma
    if i == 1
        ax = onedparpanel.sig.axes;
    else
        ax = twodparpanel.sig.axes;
    end
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    counts = nan(4,poppanel.nbins);
    for n = 1:numel(cardirs)
        L = selectedL & angL(:,n);
        counts(n,:) = histcounts(1./poppanel.(whichnd).params(L,2),ndparpanel.sig.binedges);
    end
    h = bar(ax,ndparpanel.sig.bincenters,counts','stacked');
    for n = 1:4
        set(h(n),'FaceColor',conpanel.carcols(n,:))
    end
    xlim(ax,[0 max(ndparpanel.sig.binedges)]);
    ylim(ax,[0 max(cat(2,ndparpanel.sig.counts.oneD,ndparpanel.sig.counts.twoD))]);
    
    % Exponent
    if i == 1
        ax = onedparpanel.exp.axes;
    else
        ax = twodparpanel.exp.axes;
    end
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    counts = nan(4,poppanel.nbins);
    for n = 1:numel(cardirs)
        L = selectedL & angL(:,n);
        counts(n,:) = histcounts(poppanel.(whichnd).params(L,5),ndparpanel.exp.binedges);
    end
    h = bar(ax,ndparpanel.exp.bincenters,counts','stacked');
    for n = 1:4
        set(h(n),'FaceColor',conpanel.carcols(n,:))
    end
    xlim(ax,[0 max(ndparpanel.exp.binedges)]);
    ylim(ax,[0 max(cat(2,ndparpanel.exp.counts.oneD,ndparpanel.exp.counts.twoD))]);
    
    % kappa
    if i == 1
        ax = onedparpanel.kappa.axes;
    else
        ax = twodparpanel.kappa.axes;
    end
    cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    counts = nan(4,poppanel.nbins);
    for n = 1:numel(cardirs)
        L = selectedL & angL(:,n);
        counts(n,:) = histcounts(poppanel.(whichnd).params(L,end),ndparpanel.kappa.binedges);
    end
    h = bar(ax,ndparpanel.kappa.bincenters,counts','stacked');
    for n = 1:4
        set(h(n),'FaceColor',conpanel.carcols(n,:))
    end
    xlim(ax,[0 max(ndparpanel.kappa.binedges)]);
    ylim(ax,[0 max(cat(2,ndparpanel.kappa.counts.oneD,ndparpanel.kappa.counts.twoD))]);
    
end

% Save figure variables
set(ParamsFig.poppanel,'userdata',poppanel)
set(ParamsFig.conpanel,'userdata',conpanel)
set(ParamsFig.ndparpanel,'userdata',ndparpanel)
set(ParamsFig.onedparpanel,'userdata',onedparpanel)
set(gcf,'userdata',ParamsFig);

end

function DispTable()
global GLMSPopData

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

% Unpack data
datatypes = GLMSPopData(1,:);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = num2cell(poppanel.diffnormLL);
data(:,4) = num2cell(([GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]')./pi*180);

% Repackage data and name columns
colname = {'Datafile' 'Sub' 'normLL' 'Tuning'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
strs = {'oneD' 'twoD'};
bkgndcols = ones(size(data,1),3);
for i = 1:2
    if i == 1
        ndL = poppanel.oneD.L;
    else
        ndL = ~poppanel.oneD.L;
    end
    for n = 1:4
        L = conpanel.fitcard.(strs{i}).angL(:,n) & ndL;
        bkgndcols(L,:) = repmat(conpanel.carcols(n,:),sum(L),1);
    end
end
bkgndcols(~poppanel.oneDL & ~poppanel.twoDL,:) = repmat([1 1 1],sum(~poppanel.oneDL & ~poppanel.twoDL),1);
bkgndcols(poppanel.excludeL,:) = repmat(.5,sum(poppanel.excludeL),3);

% Display table
conpanel.table = uitable('parent',ParamsFig.conpanel,...
    'units','normalized','pos',[.01 .6  .98 .39],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',bkgndcols,...
    'cellselectioncallback',@cellselect);

% Save user data
set(ParamsFig.conpanel,'userdata',conpanel)
set(ParamsFig.poppanel,'userdata',poppanel)
set(552,'userdata',ParamsFig)

end

function resetplots(~,~)

% Reset plots
UnpackPopulationData()
pulloutpeaks()
PlotTuningHists()
PlotParamHists()
DispTable()
plotnormLL()

end


%%% Interactive features %%%

function cellselect(~,b)
global GLMSPopData

% First, reset the plots
resetplots()

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndparpanel = get(ParamsFig.ndparpanel,'userdata');
%onedparpanel = get(ParamsFig.onedparpanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if ~isempty(b.Indices)
    conpanel.selectedidx = b.Indices(1);
    conpanel.popidx = conpanel.selectedidx+1;
else
    return
end

% Unpack variables
datatypes = GLMSPopData(1,:);
nd = str2double(GLMSPopData{conpanel.popidx,strcmp(datatypes,'nD')}(1));
GLMP = GLMSPopData{conpanel.popidx,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Subunit')};
Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;
surfparams = poppanel.params(conpanel.selectedidx,:);

%%% Plot Surface %%%
if nd == 1
    whichnd = 'oneD';
    selectedL = poppanel.oneDL;
elseif nd == 2
    whichnd = 'twoD';
    selectedL = poppanel.twoDL;
end
% colormap('cool')
% x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
% [xx,yy] = meshgrid(x,x);
% surface = ComputeNakaRushtonJPW(surfparams,[xx(:) yy(:)],poppanel.surftype);
% surface = reshape(surface,size(xx));
% axes(conpanel.axes.surf); cla; hold on; grid on;
% ticks = linspace(min(Lcc),max(Lcc),5);
% set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
% p = surf(gca,xx,yy,surface);
% set(p,'edgecolor','none')
% alpha(.3);
% contour3(xx,yy,surface,'linewidth',2);
% uniquestim = unique([Lcc Mcc],'rows');
% maxnsp = max(GLMP.subunit{sub}.meannspikes);
% for i = 1:size(uniquestim,1)
%     L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
%     mn = mean(nsp(L))/maxnsp*10;
%     h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
%     set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
% end
% xlabel('Lcc'); ylabel('Mcc'); zlabel('# of spikes'); title('Surface Fit')


% Only highlight histograms if the datafile is not excluded
if ~any(find(poppanel.excludeL) == conpanel.selectedidx)
    
    %%% Highlight bins in hist %%%
    % Upper A
    parval = poppanel.(whichnd).params(conpanel.selectedidx,1);
    idx = find(histcounts(parval,ndparpanel.upperA.binedges));
    binedges = ndparpanel.upperA.binedges(idx:idx+1);
    count = ndparpanel.upperA.counts.(whichnd)(idx);
    axes(ndparpanel.upperA.axes)
    plot([binedges(1) binedges(1) binedges(2) binedges(2) binedges(1)],...
        [0 count count 0 0],'r','LineWidth',2)
    
    % Baseline
    parval = poppanel.(whichnd).params(conpanel.selectedidx,6);
    idx = find(histcounts(parval,ndparpanel.bl.binedges));
    binedges = ndparpanel.bl.binedges(idx:idx+1);
    count = ndparpanel.bl.counts.(whichnd)(idx);
    axes(ndparpanel.bl.axes)
    plot([binedges(1) binedges(1) binedges(2) binedges(2) binedges(1)],...
        [0 count count 0 0],'r','LineWidth',2)
    
    % Sigma
    parval = 1/poppanel.(whichnd).params(conpanel.selectedidx,2);
    idx = find(histcounts(parval,ndparpanel.sig.binedges));
    binedges = ndparpanel.sig.binedges(idx:idx+1);
    count = ndparpanel.sig.counts.(whichnd)(idx);
    axes(ndparpanel.sig.axes)
    plot([binedges(1) binedges(1) binedges(2) binedges(2) binedges(1)],...
        [0 count count 0 0],'r','LineWidth',2)
    
    % Exponent
    parval = poppanel.(whichnd).params(conpanel.selectedidx,5);
    idx = find(histcounts(parval,ndparpanel.exp.binedges));
    binedges = ndparpanel.exp.binedges(idx:idx+1);
    count = ndparpanel.exp.counts.(whichnd)(idx);
    axes(ndparpanel.exp.axes)
    plot([binedges(1) binedges(1) binedges(2) binedges(2) binedges(1)],...
        [0 count count 0 0],'r','LineWidth',2)
    
    % Kappa
    parval = poppanel.(whichnd).params(conpanel.selectedidx,end);
    idx = find(histcounts(parval,ndparpanel.kappa.binedges));
    binedges = ndparpanel.kappa.binedges(idx:idx+1);
    count = ndparpanel.kappa.counts.(whichnd)(idx);
    axes(ndparpanel.kappa.axes)
    plot([binedges(1) binedges(1) binedges(2) binedges(2) binedges(1)],...
        [0 count count 0 0],'r','LineWidth',2)
    
    %%% Highlight bins in pop tuning %%%
    % Find the selected collection of cells
    ang = surfparams(end-1);
    axes(poppanel.axes.(whichnd))
    [tout,~] = rose(poppanel.tuning);
    edges = unique(tout)-pi;
    lowerTheta = edges(find(sign(ang-edges)==1,1,'last'));
    upperTheta = edges(find(sign(ang-edges)==-1,1,'first'));
    angL = poppanel.tuning >= lowerTheta & poppanel.tuning <= upperTheta;
    L = selectedL & angL;
    p = rose(poppanel.tuning(L));
    set(p,'color','r','markerfacecolor',[1 0 0])
    
end

%%% Highlight cell in table %%%
conpanel.table.BackgroundColor(conpanel.selectedidx,:) = [1 0 0];

% Save user data
set(ParamsFig.poppanel,'UserData',poppanel)
set(ParamsFig.ndparpanel,'UserData',ndparpanel)
set(ParamsFig.conpanel,'userdata',conpanel)
set(gcf,'UserData',ParamsFig)


end


%%% Call other analyses %%%

function callanalysis(~,~)
global GLMSPopData GLMP

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');

% If no cell selected, return prompt. If selected, call analysis.
if isempty(conpanel.selectedidx)
    disp('Must select a cell first')
else
    datatype = GLMSPopData(1,:);
    GLMP = GLMSPopData{conpanel.popidx,strcmp(datatype,'GLMP')};
    sub = GLMSPopData{conpanel.popidx,strcmp(datatype,'Subunit')};
    GLMSGUI_Surface([],[],[],sub)
end

% Save user data
set(ParamsFig.conpanel,'userdata',conpanel)
set(552,'userdata',ParamsFig)

cellselect()

end

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
DN = GLMSPopData{conpanel.popidx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{conpanel.popidx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end


%%% Setup %%%

function setUpFig()

figure(552); clf;
set(gcf,'units','normalized','pos',[.05 .1 .9 .8],'NumberTitle','off','Name','Population Parameters')

% Set up panels
ParamsFig.conpanel = uipanel('units','normalized','pos',[.01 .01 .19 .98],'title','Control Panel');
ParamsFig.poppanel = uipanel('units','normalized','pos',[.21 .01 .18 .98],'title','Population Tuning');
ParamsFig.ndparpanel = uipanel('units','normalized','pos',[.4 .01 .19 .98],'title','1D vs 2D');
ParamsFig.onedparpanel = uipanel('units','normalized','pos',[.6 .01 .19 .98],'title','1D Params');
ParamsFig.twodparpanel = uipanel('units','normalized','pos',[.8 .01 .19 .98],'title','2D Params');

% Set up population panel
poppanel.axes.oneD = axes('parent',ParamsFig.poppanel,'unit','normalized',...
    'pos',[.1 .575 .8 .4]);
poppanel.axes.twoD = axes('parent',ParamsFig.poppanel,'units','normalized',...
    'pos',[.1 .15 .8 .4]);
strs = {'All 1D Surfs' 'Unichromatic' 'Bichromatic'};
poppanel.uicontrols.which1Dsurf = uicontrol(ParamsFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .575 .8 .02],...
    'string',strs,'value',1,'callback',@resetplots);
strs = {'All 2D Surfs' 'Unichromatic hyperbola' 'Bichromatic hyperbola' 'Unichromatic elipse' 'Bichromatic elipse'};
poppanel.uicontrols.which2Dsurf = uicontrol(ParamsFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .125 .8 .02],...
    'string',strs,'value',1,'callback',@resetplots);
strs = {'Both Monkeys' 'Nut Data' 'Maui Data'};
poppanel.uicontrols.whichMonk = uicontrol('parent',ParamsFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .05 .8 .05],...
    'string',strs,'callback',@resetplots);
poppanel.nbins = 20;
poppanel.oneDL = [];
poppanel.twoDL = [];

% Set up Control Panel
conpanel.table = uitable('parent',ParamsFig.conpanel,...
    'units','normalized','pos',[.01 .6  .98 .39]);
conpanel.axes.hist = axes('parent',ParamsFig.conpanel,'units','normalized',...
    'pos',[.15 .3 .75 .25]); box on;
conpanel.uicontrols.overview = uicontrol('parent',ParamsFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .01 .4 .05],'callback',@LoadOverview);
conpanel.uicontrols.analindiv = uicontrol('parent',ParamsFig.conpanel,...
    'style','pushbutton','string','Surface Analysis','units','normalized',...
    'pos',[.525 .01 .4 .05],'callback',@callanalysis);
conpanel.uicontrols.wedgesize = uicontrol(ParamsFig.conpanel,'style','edit',...
    'units','normalized','pos',[.6 .07 .3 .04],'string',45,...
    'callback',@resetplots);
conpanel.uicontrols.wedgesizelabel = uicontrol(ParamsFig.conpanel,'style','text',...
    'units','normalized','pos',[.6 .11 .3 .04],'string','Wedge Size');
conpanel.uicontrols.ndthresh = uicontrol(ParamsFig.conpanel,'style','edit',...
    'units','normalized','pos',[.6 .2 .3 .04],'string',.075,...
    'callback',@resetplots);
conpanel.uicontrols.ndthreshlabel = uicontrol(ParamsFig.conpanel,'style','text',...
    'units','normalized','pos',[.1 .2 .45 .04],'string','Norm LL 1D/2D threshold:');
conpanel.selectedidx = [];
conpanel.popidx = [];

% Some hard coded values
conpanel.oneDcol = [1 .5 0];
conpanel.twoDcol = [0 0 1];
conpanel.cardirs = [pi/4 -3*pi/4 -pi/4 3*pi/4];
conpanel.carcols = [.8 .8 .5; .1 .1 .5; 1 .5 .5; .5 1 .5];

% Set up par panel
ndparpanel = [];

% Save figure variables
set(ParamsFig.conpanel,'userdata',conpanel)
set(ParamsFig.poppanel,'userdata',poppanel)
set(ParamsFig.ndparpanel,'userdata',ndparpanel)
set(gcf,'userdata',ParamsFig)

end

function SetUpHists()

% Load figure variables
ParamsFig = get(gcf,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndparpanel = get(ParamsFig.ndparpanel,'userdata');
%onedparpanel = get(ParamsFig.onedparpanel,'userdata');
%twodparpanel = get(ParamsFig.twodparpanel,'userdata');

% Set up Parameter Space
% Upper asmyptote
ndparpanel.upperA.axes = axes('parent',ParamsFig.ndparpanel,'unit','normalized',...
    'pos',[.1 .825 .8 .15],'box','on','title','Upper Asymptote');
ndparpanel.upperA.labels = uicontrol('style','text','parent',ParamsFig.ndparpanel,...
    'units','normalized','position',[.1 .975 .8 .02],'string','Upper Asymptote',...
    'fontsize',12,'horizontalalignment','center');
onedparpanel.upperA.axes = axes('parent',ParamsFig.onedparpanel,'unit','normalized',...
    'pos',[.1 .825 .8 .15],'box','on','title','Upper Asymptote');
onedparpanel.upperA.labels = uicontrol('style','text','parent',ParamsFig.onedparpanel,...
    'units','normalized','position',[.1 .975 .8 .02],'string','Upper Asymptote',...
    'fontsize',12,'horizontalalignment','center');
twodparpanel.upperA.axes = axes('parent',ParamsFig.twodparpanel,'unit','normalized',...
    'pos',[.1 .825 .8 .15],'box','on','title','Upper Asymptote');
twodparpanel.upperA.labels = uicontrol('style','text','parent',ParamsFig.twodparpanel,...
    'units','normalized','position',[.1 .975 .8 .02],'string','Upper Asymptote',...
    'fontsize',12,'horizontalalignment','center');

% Basline
ndparpanel.bl.axes = axes('parent',ParamsFig.ndparpanel,'units','normalized',...
    'pos',[.1 .625 .8 .15],'box','on');
ndparpanel.bl.labels = uicontrol('style','text','parent',ParamsFig.ndparpanel,...
    'units','normalized','position',[.1 .775 .8 .02],'string','Baseline',...
    'fontsize',12,'horizontalalignment','center');
onedparpanel.bl.axes = axes('parent',ParamsFig.onedparpanel,'units','normalized',...
    'pos',[.1 .625 .8 .15],'box','on');
onedparpanel.bl.labels = uicontrol('style','text','parent',ParamsFig.onedparpanel,...
    'units','normalized','position',[.1 .775 .8 .02],'string','Baseline',...
    'fontsize',12,'horizontalalignment','center');
twodparpanel.bl.axes = axes('parent',ParamsFig.twodparpanel,'units','normalized',...
    'pos',[.1 .625 .8 .15],'box','on');
twodparpanel.bl.labels = uicontrol('style','text','parent',ParamsFig.twodparpanel,...
    'units','normalized','position',[.1 .775 .8 .02],'string','Baseline',...
    'fontsize',12,'horizontalalignment','center');

% Sigmas
ndparpanel.sig.axes = axes('parent',ParamsFig.ndparpanel,'unit','normalized',...
    'pos',[.1 .425 .8 .15],'box','on');
ndparpanel.sig.labels = uicontrol('style','text','parent',ParamsFig.ndparpanel,...
    'units','normalized','position',[.1 .575 .8 .02],'string','Sigma',...
    'fontsize',12,'horizontalalignment','center');
onedparpanel.sig.axes = axes('parent',ParamsFig.onedparpanel,'unit','normalized',...
    'pos',[.1 .425 .8 .15],'box','on');
onedparpanel.sig.labels = uicontrol('style','text','parent',ParamsFig.onedparpanel,...
    'units','normalized','position',[.1 .575 .8 .02],'string','Sigma',...
    'fontsize',12,'horizontalalignment','center');
twodparpanel.sig.axes = axes('parent',ParamsFig.twodparpanel,'unit','normalized',...
    'pos',[.1 .425 .8 .15],'box','on');
twodparpanel.sig.labels = uicontrol('style','text','parent',ParamsFig.twodparpanel,...
    'units','normalized','position',[.1 .575 .8 .02],'string','Sigma',...
    'fontsize',12,'horizontalalignment','center');

% Exponent
ndparpanel.exp.axes = axes('parent',ParamsFig.ndparpanel,'unit','normalized',...
    'pos',[.1 .225 .8 .15],'box','on');
ndparpanel.exp.labels = uicontrol('style','text','parent',ParamsFig.ndparpanel,...
    'units','normalized','position',[.1 .375 .8 .02],'string','Exponent',...
    'fontsize',12,'horizontalalignment','center');
onedparpanel.exp.axes = axes('parent',ParamsFig.onedparpanel,'unit','normalized',...
    'pos',[.1 .225 .8 .15],'box','on');
onedparpanel.exp.labels = uicontrol('style','text','parent',ParamsFig.onedparpanel,...
    'units','normalized','position',[.1 .375 .8 .02],'string','Exponent',...
    'fontsize',12,'horizontalalignment','center');
twodparpanel.exp.axes = axes('parent',ParamsFig.twodparpanel,'unit','normalized',...
    'pos',[.1 .225 .8 .15],'box','on');
twodparpanel.exp.labels = uicontrol('style','text','parent',ParamsFig.twodparpanel,...
    'units','normalized','position',[.1 .375 .8 .02],'string','Exponent',...
    'fontsize',12,'horizontalalignment','center');

% Kappa
ndparpanel.kappa.axes = axes('parent',ParamsFig.ndparpanel,'unit','normalized',...
    'pos',[.1 .025 .8 .15],'box','on');
ndparpanel.kappa.labels = uicontrol('style','text','parent',ParamsFig.ndparpanel,...
    'units','normalized','position',[.1 .175 .8 .02],'string','Kappa',...
    'fontsize',12,'horizontalalignment','center');
onedparpanel.kappa.axes = axes('parent',ParamsFig.onedparpanel,'unit','normalized',...
    'pos',[.1 .025 .8 .15],'box','on');
onedparpanel.kappa.labels = uicontrol('style','text','parent',ParamsFig.onedparpanel,...
    'units','normalized','position',[.1 .175 .8 .02],'string','Kappa',...
    'fontsize',12,'horizontalalignment','center');
twodparpanel.kappa.axes = axes('parent',ParamsFig.twodparpanel,'unit','normalized',...
    'pos',[.1 .025 .8 .15],'box','on');
twodparpanel.kappa.labels = uicontrol('style','text','parent',ParamsFig.twodparpanel,...
    'units','normalized','position',[.1 .175 .8 .02],'string','Kappa',...
    'fontsize',12,'horizontalalignment','center');

% Save Data
set(ParamsFig.poppanel,'UserData',poppanel)
set(ParamsFig.ndparpanel,'UserData',ndparpanel)
set(ParamsFig.onedparpanel,'UserData',onedparpanel)
set(ParamsFig.twodparpanel,'userdata',twodparpanel)
set(gcf,'UserData',ParamsFig)

end

function UnpackPopulationData()
global GLMSPopData

% Load figure variables
ParamsFig = get(gcf,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndparpanel = get(ParamsFig.ndparpanel,'userdata');
thresh = str2double(get(conpanel.uicontrols.ndthresh,'String'));

% Pull out some tuning data from population
desinfo = GLMSPopData(1,:);

% Datafiles
poppanel.filenames = GLMSPopData(2:end,strcmp(desinfo,'Datafile'));

% Which monkey
monks = cat(1,poppanel.filenames{:});
poppanel.monk.NutL = monks(:,1) == 'N';
poppanel.monk.MauiL = monks(:,1) == 'M';
conpanel.monkL = true(size(poppanel.monk.NutL));

% Surf parameters
surfparams = [GLMSPopData{2:end,strcmp(desinfo,'Surface Parameters')}]';
poppanel.surftype = surfparams(1).oneD.surftype;
fitcomps = [surfparams.fitcomps]';
oneD = [fitcomps.oneD]';
parvals = cat(1,oneD.parvals);
poppanel.oneD.unichrom.parvals = parvals(1:2:end,:);
poppanel.oneD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [oneD.nLLs]';
poppanel.oneD.unichrom.LLs = -nLLs(:,1);
poppanel.oneD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,oneD.normLL);
poppanel.oneD.unichrom.normLLs = normLLs(1:2:end);
poppanel.oneD.bichrom.normLLs = normLLs(2:2:end);

twoD = [fitcomps.twoD]';
parvals = cat(1,twoD.parvals);
poppanel.twoD.unichrom.parvals = parvals(1:2:end,:);
poppanel.twoD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [twoD.nLLs]';
poppanel.twoD.unichrom.LLs = -nLLs(:,1);
poppanel.twoD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,twoD.normLL);
poppanel.twoD.unichrom.normLLs = normLLs(1:2:end);
poppanel.twoD.bichrom.normLLs = normLLs(2:2:end);

% Distinguish 1D from 2D using bichromatic fits
poppanel.oneD.L = (poppanel.twoD.bichrom.normLLs - poppanel.oneD.bichrom.normLLs) < thresh;

% Distinguish unichrom from bichrom
poppanel.oneD.unichrom.L = (poppanel.oneD.bichrom.normLLs - poppanel.oneD.unichrom.normLLs) < thresh;
poppanel.oneD.bichrom.L = (poppanel.oneD.bichrom.normLLs - poppanel.oneD.unichrom.normLLs) >= thresh;
poppanel.twoD.unichrom.L = (poppanel.twoD.bichrom.normLLs - poppanel.twoD.unichrom.normLLs) < thresh;
poppanel.twoD.bichrom.L = (poppanel.twoD.bichrom.normLLs - poppanel.twoD.unichrom.normLLs) >= thresh;

% Unichrom and bichrom parvals
poppanel.oneD.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.oneD.params(poppanel.oneD.unichrom.L,:) = poppanel.oneD.unichrom.parvals(poppanel.oneD.unichrom.L,:);
poppanel.oneD.params(poppanel.oneD.bichrom.L,:) = poppanel.oneD.bichrom.parvals(poppanel.oneD.bichrom.L,:);
poppanel.twoD.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.twoD.params(poppanel.twoD.unichrom.L,:) = poppanel.twoD.unichrom.parvals(poppanel.twoD.unichrom.L,:);
poppanel.twoD.params(poppanel.twoD.bichrom.L,:) = poppanel.twoD.bichrom.parvals(poppanel.twoD.bichrom.L,:);

% All parvals
poppanel.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.params(poppanel.oneD.L,:) = poppanel.oneD.params(poppanel.oneD.L,:);
poppanel.params(~poppanel.oneD.L,:) = poppanel.twoD.params(~poppanel.oneD.L,:);
   
% Sort by surface type
orax = poppanel.twoD.params(:,4);
poppanel.twoD.unichrom.hypL = orax < 0 & poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.unichrom.eliL = orax > 0 & poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.bichrom.hypL = orax < 0 & ~poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.bichrom.eliL = orax > 0 & ~poppanel.twoD.unichrom.L & ~poppanel.oneD.L;

% Tuning
poppanel.tuning = nan(size(poppanel.oneD.L));
poppanel.tuning(poppanel.oneD.L) = poppanel.oneD.params(poppanel.oneD.L,end-1);
poppanel.tuning(~poppanel.oneD.L) = poppanel.twoD.params(~poppanel.oneD.L,end-1);
poppanel.diffnormLL = poppanel.twoD.bichrom.normLLs - poppanel.oneD.bichrom.normLLs;


% Save figure variables
set(ParamsFig.poppanel,'userdata',poppanel)
set(ParamsFig.conpanel,'userdata',conpanel)
set(ParamsFig.ndparpanel,'userdata',ndparpanel)
set(gcf,'userdata',ParamsFig);

exclusioncriteria();

end

function exclusioncriteria()
global GLMSPopData

% Load figure variables
ParamsFig = get(gcf,'userdata');
%conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

% parmeter bounds
modmin = 10; % in sp/s
kappamax = 4;
%confmax = 360;

% Load data
datatypes = GLMSPopData(1,:);

% Modulation
modvals = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    modvals(n-1) = max(GLMP.subunit{sub}.meannspikes)...
        ./ mean(GLMP.subunit{sub}.stimDur) - mean(GLMP.subunit{sub}.blfr);
end
poppanel.modL = modvals >= modmin;
disp([num2str(sum(~poppanel.modL)) ' datafiles do not meet modulation criteria.'])

% Kappa
kappavals = nan(size(GLMSPopData,1)-1,1);
kappavals(poppanel.oneD.L) = poppanel.oneD.params(poppanel.oneD.L,end);
kappavals(~poppanel.oneD.L) = poppanel.twoD.params(~poppanel.oneD.L,end);
poppanel.kappaL = kappavals <= kappamax;
disp([num2str(sum(~poppanel.kappaL)) ' datafiles do not meet variance criteria.'])


% Radially symmmetric datasets
RSL = zeros(size(modvals)); 
RSL(1:28) = 1; % hard coded, last RS dataset
poppanel.RSL = RSL;
disp([num2str(sum(poppanel.RSL)) ' datafiles are RS.'])


% Which don't meet criteria?
poppanel.excludeL = modvals < modmin | kappavals > kappamax | RSL;% | poppanel.conf > confmax;
% disp([num2str(sum(poppanel.excludeL)) ' datafiles excluded based on current criteria.'])

% For selecting one subunit per datafile
filenames = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[~,fnidx,ufnidx] = unique(filenames); 
repsL = ones(size(filenames,1),1); 
repsL(fnidx) = 0; %for 1 sub/cell
repeatsidx = find(repsL);
poppanel.uniqL = ones(size(filenames,1),1);
for n = 1:numel(repeatsidx)
    idx = find(ufnidx == ufnidx(repeatsidx(n)));
    glmp1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'GLMP')};
    sub1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'Subunit')};
    glmp2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'GLMP')};
    sub2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'Subunit')};
    maxnsp1 = max(glmp1.subunit{sub1}.meannspikes);
    maxnsp2 = max(glmp2.subunit{sub2}.meannspikes);
    if any(poppanel.excludeL(idx))
        poppanel.uniqL(idx) = ~poppanel.excludeL(idx);
    elseif maxnsp1 > maxnsp2
        poppanel.uniqL(idx(2)) = 0;
    else
        poppanel.uniqL(idx(1)) = 0;
    end
end
disp([num2str(sum(~poppanel.uniqL)) ' secondary subunits excluded from analysis.'])

poppanel.excludeL = poppanel.excludeL | ~poppanel.uniqL;
disp([num2str(sum(poppanel.excludeL)) '/' num2str(numel(poppanel.excludeL)) ' total datafiles excluded from analysis.'])

% Save user data
set(ParamsFig.poppanel,'userdata',poppanel)
set(gcf,'userdata',ParamsFig)


end

function plotnormLL()

% Load figure variables
ParamsFig = get(gcf,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

thresh = str2double(conpanel.uicontrols.ndthresh.String);
L = ~poppanel.excludeL;
diffnormLL = poppanel.twoD.bichrom.normLLs(L) - poppanel.oneD.bichrom.normLLs(L);

% populate normLL hist
axes(conpanel.axes.hist); cla; hold on;
maxval = max(diffnormLL);
%h = histogram(poppanel.diffnormLL,(0:.005:maxval),'edgecolor','w','facecolor','k');
h = histogram(diffnormLL,linspace(0,maxval,50),'edgecolor','w','facecolor','k');
ylim([0 max(h.Values)*1.1])
xlim([0 (rndofferr(maxval,2))])
plot([thresh thresh],[0 max(h.Values)*1.1],'r--')


% Save user data
set(ParamsFig.poppanel,'userdata',poppanel)
set(ParamsFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',ParamsFig)


end

function pulloutpeaks()

% Load figure variables
ParamsFig = get(gcf,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
pmdeg = str2double(get(conpanel.uicontrols.wedgesize,'string'))/2;
pmrad = pmdeg/180*pi;


% Find peaks of PD distributions
angs1d = poppanel.tuning(poppanel.oneD.L & conpanel.monkL & ~poppanel.excludeL);
angs2d = poppanel.tuning(~poppanel.oneD.L & conpanel.monkL & ~poppanel.excludeL);

% Create histogram and smooth
nbins = 361;
bins = linspace(-pi,pi,nbins);
binsize = mean(diff(bins));
gaussSize = 45.001*pi/180; %in rad (ugh so stupid, at 45 the -L-M is off, but at this SLIGHT amount more, its right on)
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
PSTH1d = hist(angs1d,bins);
PSTH2d = hist(angs2d,bins);
smoothPSTH = conv(PSTH1d,gaussfilter,'same');
smoothPSTH1d = smoothPSTH./max(smoothPSTH) .* max(PSTH1d); % normalizing by max height of PSTH for plotting purposes
smoothPSTH = conv(PSTH2d,gaussfilter,'same');
smoothPSTH2d = smoothPSTH./max(smoothPSTH) .* max(PSTH2d); % normalizing by max height of PSTH for plotting purposes

A1 = spline(bins,smoothPSTH1d);
A2 = spline(bins,smoothPSTH2d);
M = diag(3:-1:1,1);

%first derivative
D1 = A1;
D1.coefs = D1.coefs*M;
D2 = A2;
D2.coefs = D2.coefs*M;
conpanel.fitcard.oneD.means = nan(4,1);
conpanel.fitcard.oneD.angL = nan(numel(poppanel.tuning),4);
conpanel.fitcard.twoD.means = nan(4,1);
conpanel.fitcard.twoD.angL = nan(numel(poppanel.tuning),4);
for n = 1:4
    conpanel.fitcard.oneD.means(n) = fzero(@(bins) ppval(D1,bins),conpanel.cardirs(n));
    conpanel.fitcard.twoD.means(n) = fzero(@(bins) ppval(D2,bins),conpanel.cardirs(n));
    conpanel.fitcard.oneD.angL(:,n) = poppanel.tuning > conpanel.fitcard.oneD.means(n)-pmrad...
        & poppanel.tuning < conpanel.fitcard.oneD.means(n)+pmrad;
    conpanel.fitcard.twoD.angL(:,n) = poppanel.tuning > conpanel.fitcard.twoD.means(n)-pmrad...
        & poppanel.tuning < conpanel.fitcard.twoD.means(n)+pmrad;
end

% Save user data
set(ParamsFig.poppanel,'userdata',poppanel)
set(ParamsFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',ParamsFig)

end
