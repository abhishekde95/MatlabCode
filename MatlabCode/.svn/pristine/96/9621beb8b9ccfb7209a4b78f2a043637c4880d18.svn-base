function GLMSPopGUI_BLSpikeStats(~,~)

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

figure(16); clf;
set(gcf,'pos',[50 75 1200 750],'numbertitle','off','name','Population Nonstationarities');

% Set up panels
BlSpFig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .4 .49]);
BlSpFig.popanel = uipanel(gcf,'units','normalized','pos',[.01 .51 .4 .48]);
BlSpFig.cellpanel = uipanel(gcf,'units','normalized','pos',[.42 .01 .57 .98]);

% Set up conpanel
conpanel.uicontrols.reanalyzeall = uicontrol('parent',BlSpFig.conpanel,...
    'style','pushbutton','string','Reanalyze All Cells','units','normalized',...
    'pos',[.05 .01 .425 .08],'Foregroundcolor',[1 0 0],'callback',@ReanalAll);
conpanel.uicontrols.analindiv = uicontrol('parent',BlSpFig.conpanel,...
    'style','pushbutton','string','PSTH Analysis','units','normalized',...
    'pos',[.525 .01 .425 .08],'callback',@callanalysis);
conpanel.uicontrols.overview = uicontrol('parent',BlSpFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .1 .425 .08],'callback',@LoadOverview);
conpanel.uicontrols.reanalyze = uicontrol('parent',BlSpFig.conpanel,...
    'style','pushbutton','string','Reanalyze Cell','units','normalized',...
    'pos',[.525 .1 .425 .08],'callback',@Reanal);
conpanel.selectedidx = [];
conpanel.axlims = [-.2 .6];
conpanel.binsize = .01; % in seconds

% Set up popanel
popanel.axes.pvals = axes('parent',BlSpFig.popanel,'units','normalized',...
    'pos',[.1 .1 .8 .8]); box on;

% Set up cellpanel
cellpanel.axes.rasters = axes('parent',BlSpFig.cellpanel,'units','normalized',...
    'pos',[.4 .4 .55 .55],'box','on');
cellpanel.axes.bubbleplot = axes('parent',BlSpFig.cellpanel,'units','normalized',...
    'pos',[.1 .1 .2 .2],'box','on'); axis square;
cellpanel.axes.bldist = axes('parent',BlSpFig.cellpanel,'units','normalized',...
    'pos',[.1 .4 .2 .55],'box','on');
cellpanel.axes.PSTH = axes('parent',BlSpFig.cellpanel,'units','normalized',...
    'pos',[.4 .1 .55 .2],'box','on');

% Save user data
set(BlSpFig.conpanel,'userdata',conpanel)
set(BlSpFig.popanel,'userdata',popanel)
set(BlSpFig.cellpanel,'userdata',cellpanel)
set(gcf,'userdata',BlSpFig)

% Unpack pop data
UnpackPopulationData()

% Display table
DispTable()

% Show stats
ShowPopData()

end

function cellselect(~,b)
global GLMSPopData

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');
popanel = get(BlSpFig.popanel,'userdata');
cellpanel = get(BlSpFig.cellpanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if ~isempty(b.Indices)
    conpanel.selectedidx = b.Indices(1);
else
    return
end
set(conpanel.table,'CellSelectionCallback',[])
popidx = conpanel.selectedidx + 1;

% Unpack variables
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{popidx,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{popidx,strcmp(datatypes,'Subunit')};

% Set up for rasters by trial
nrows = numel(GLMP.subunit{sub}.Lcc);
rowcoords = linspace(1,0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% Plot PSTH
bins = conpanel.axlims(1):conpanel.binsize:conpanel.axlims(2);
PSTH = histc(cat(1,GLMP.subunit{sub}.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.subunit{sub}.Lcc)) ./ conpanel.binsize;
axes(cellpanel.axes.PSTH); cla; hold on;
bar(bins,PSTH,'facecolor',[.5 0 .9],'edgecolor','k');

set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

% Plot rasters
axes(cellpanel.axes.rasters); cla; hold on;
title('Rasters by Trial')
xlim(conpanel.axlims)
ylim([0 1])
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 1],'m-')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 1],'m-')
for n = 1:nrows
    tpts = GLMP.subunit{sub}.normspiketimes{n};
    if ~isempty(tpts)
        plot(repmat([tpts tpts],numel(tpts),1),repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts,1)),'k-')
    end
end
set(gca,'YTick',[0 .25 .5 .75 1],'YTickLabel',[nrows round(3*nrows/4) round(nrows/2) round(nrows/4) 1])
xlabel('Normalized Spike Time')
ylabel('Trial Number')

% Bubble plot
axes(cellpanel.axes.bubbleplot); cla; hold on; grid on; axis square;
maxnsp = max(GLMP.subunit{sub}.meannspikes);
Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;
uniquestim = unique([Lcc Mcc],'rows');
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('Bubble Plot')

% Plot rank order
axes(cellpanel.axes.bldist); cla; box on; hold on;
trs = (1:numel(GLMP.subunit{sub}.blfr))';
blfrs = GLMP.subunit{sub}.blfr;
%[rho,pval] = corr(trs,blfrs,'type','Spearman');
plot(blfrs,trs,'ko')
set(gca,'YDir','reverse')
ylabel('Firing Rate (sp/s)')
ylim([1 max(trs)])
xlim([0 max(blfrs)+1])

% Highlight bin in pop
% bins = linspace(0,1,20);
% axes(popanel.axes.pvals); cla; hold on;
% [counts] = hist(conpanel.pvals,bins);
% tempvect = zeros(size(counts));
% tempvect(ind(conpanel.selectedidx)) = counts(ind(conpanel.selectedidx));
% bar(bins,counts)
% bar(bins,tempvect,'r')

set(BlSpFig.cellpanel,'title',[GLMP.datafile ' sub # ' num2str(sub)]);

% Save user data
set(conpanel.table,'CellSelectionCallback',@cellselect)
set(BlSpFig.conpanel,'userdata',conpanel)
set(BlSpFig.cellpanel,'userdata',cellpanel)
set(BlSpFig.popanel,'userdata',popanel)
set(gcf,'userdata',BlSpFig)

end

function ShowPopData()

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');
popanel = get(BlSpFig.popanel,'userdata');

% Plot pop pvals (and rhos?)
axes(popanel.axes.pvals); cla; hold on;
%[popanel.counts,popanel.cents] = hist(conpanel.pvals,20);
hist(conpanel.pvals,20);
%plot(conpanel.pvals,conpanel.rhos,'ko')

% Save user data
set(BlSpFig.popanel,'userdata',popanel)
set(16,'userdata',BlSpFig)

end

function ReanalAll(~,~,~)
global GLMSPopData

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');

% Rotate through selected index
for n = 1:(size(GLMSPopData,1)-1)
    conpanel.selectedidx = n;
    set(BlSpFig.conpanel,'userdata',conpanel)
    set(gcf,'userdata',BlSpFig)
    Reanal()
end

disp('Finished reanalyzing baseline rank correlation.')

end

function Reanal(~,~,~)
global GLMSPopData GLMP

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');

% Which dataset
if isempty(conpanel.selectedidx);
    disp('Must Select Cell First')
    return
end

% Reload nex file, unpack, analyze, and save
datatypes = GLMSPopData(1,:);
if ismac
    library1 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
    library2 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library1 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
    library2 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
datafile = conpanel.table.Data(conpanel.selectedidx,1);
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
disp(['Analyzing datafile No. ' num2str(conpanel.selectedidx) ': ' datafile{1} ' sub # ' num2str(sub)])

% Unpack data and begin surface analysis
popidx =  conpanel.selectedidx+1;
GLMP = GLMSPopData{popidx,strcmp(datatypes,'GLMP')};
spikestats = GLMSPopData{popidx,strcmp(datatypes,'GLMPSpikeStats')};

% Set up for BL psth
trs = (1:numel(GLMP.subunit{sub}.blfr))';
blfrs = GLMP.subunit{sub}.blfr;
[spikestats.BLstats.spearrank.rho,spikestats.BLstats.spearrank.pval] = corr(trs,blfrs,'type','Spearman');

% Organize parameters into population structure and save
GLMSPopData{popidx,strcmp(datatypes,'GLMPSpikeStats')} = spikestats;

% Save user data
save([library2 'GLMSPopData'],'GLMSPopData')
set(BlSpFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',BlSpFig)

% Now that new data is in the pop structure, reload whole analysis
UnpackPopulationData()
DispTable()
b.Indices = conpanel.selectedidx;
cellselect([],b)

disp(['Datafile No. ' num2str(conpanel.selectedidx) ': ' datafile{1} ' sub # ' num2str(sub) ' Reanalyzed.'])

end

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
if isempty(idx)
    disp('Must select a datafile before overview can be shown.')
    return
end
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

function callanalysis(~,~)
global GLMSPopData GLMP

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');

% If no cell selected, return prompt. If selected, call analysis.
if isempty(conpanel.selectedidx)
    disp('Must select a cell first')
else
    datatype = GLMSPopData(1,:);
    GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatype,'GLMP')};
    sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatype,'Subunit')};
    GLMSGUI_PSTH([],[],[],sub)
end

% Save user data
set(BlSpFig.conpanel,'userdata',conpanel)
set(16,'userdata',BlSpFig)

b.Indices = conpanel.selectedidx; 
cellselect([],b)

end

function UnpackPopulationData()
global GLMSPopData

% Load figure variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');
popanel = get(BlSpFig.popanel,'userdata');
datatypes = GLMSPopData(1,:);

% nd params
nd = GLMSPopData(2:end,strcmp(datatypes,'nD'));
conpanel.oneD.L = strcmp(nd,'1D');
conpanel.twoD.L = strcmp(nd,'2D');

% BL Rank Correlation
temp = [GLMSPopData{2:end,strcmp(datatypes,'GLMPSpikeStats')}];
temp2 = [temp.BLstats];
temp3 = [temp2.spearrank];
conpanel.pvals = [temp3.pval]';
conpanel.rhos = [temp3.rho]';

% tuning 
conpanel.tuning = ([GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]')./180*pi;

% Monkey
filenames = cat(1,GLMSPopData{2:end,1});
monk = filenames(:,1);
conpanel.monk.M = monk == 'M';
conpanel.monk.N = monk == 'N';

% Save user data
set(BlSpFig.conpanel,'userdata',conpanel)
set(BlSpFig.popanel,'userdata',popanel)
set(gcf,'userdata',BlSpFig)

end

function DispTable()
global GLMSPopData

% Load figure and pop variables
BlSpFig = get(gcf,'userdata');
conpanel = get(BlSpFig.conpanel,'userdata');
popanel = get(BlSpFig.popanel,'userdata');

% Unpack data
datatypes = GLMSPopData(1,:);
data = cell(size(GLMSPopData,1)-1,4);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = num2cell(conpanel.pvals);
data(:,4) = num2cell(conpanel.rhos);
data(:,5) = GLMSPopData(2:end,strcmp(datatypes,'nD'));
data(:,6) = GLMSPopData(2:end,strcmp(datatypes,'Tuning'));

% Repackage data and name columns
colname = {'Datafile' 'Sub' 'P-Vals' 'Rhos' 'nD' 'Tuning'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table
conpanel.table = uitable('parent',BlSpFig.conpanel,...
    'units','normalized','pos',[.01 .2  .98 .79],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',[1 1 1],...
    'CellSelectionCallback',@cellselect);

% Save user data
set(BlSpFig.conpanel,'userdata',conpanel)
set(BlSpFig.popanel,'userdata',popanel)
set(gcf,'userdata',BlSpFig)

end



