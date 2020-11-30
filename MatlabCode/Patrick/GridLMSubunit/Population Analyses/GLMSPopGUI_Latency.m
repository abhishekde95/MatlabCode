function GLMSPopGUI_Latency(~,~)
%global GLMSPopData
% This function compares the latency measurements of each dataset

figure(51015); clf; 
set(gcf,'units','normalized','pos',[.1 .2 .75 .5],...
    'name','GLMS Population Analysis: Latency','numbertitle','off')
PopLatFig.conpanel = uipanel('units','normalized','pos',[.01 .01 .35 .98]);
PopLatFig.popanel = uipanel('units','normalized','pos',[.36 .01 .34 .98]);
PopLatFig.dispanel = uipanel('units','normalized','pos',[.71 .01 .28 .98]);

% Set up popanel
popanel.axes.oneD = axes('parent',PopLatFig.popanel,'units','normalized',...
    'pos',[.15 .6 .8 .35],'box','on','ButtonDownFcn',@highlightcells);
title('1D Subunits')
popanel.axes.twoD = axes('parent',PopLatFig.popanel,'units','normalized',...
    'pos',[.15 .1 .8 .35],'box','on','ButtonDownFcn',@highlightcells);
title('2D Subunits')

% Set up dispanel
dispanel.axes.FSP = axes('parent',PopLatFig.dispanel,'units','normalized',...
    'pos',[.15 .1 .75 .35],'box','on');
title('First Spike')
dispanel.axes.GTF = axes('parent',PopLatFig.dispanel,'units','normalized',...
    'pos',[.15 .6 .75 .35],'box','on');
title('Gaussian Temporal Window')

% Set up conpanel
conpanel.table = uitable('parent',PopLatFig.conpanel,...
    'units','normalized','pos',[.01 .3 .98 .69],...
    'BackgroundColor',[1 1 1],...
    'cellselectioncallback',@cellselect);
conpanel.uicontrols.reanalyze = uicontrol('parent',PopLatFig.conpanel,...
    'style','pushbutton','string','Reanalyze All Cells','units','normalized',...
    'pos',[.05 .05 .25 .08],'callback',@Reanal);
conpanel.uicontrols.analindiv = uicontrol('parent',PopLatFig.conpanel,...
    'style','pushbutton','string','Latency Analysis','units','normalized',...
    'pos',[.35 .05 .25 .08],'callback',@callanalysis);
str = {'First Spike' 'Gaussian Temporal Form'};
conpanel.uicontrols.latanalselector = uicontrol(PopLatFig.conpanel,...
    'style','popupmenu','string',str,'units','normalized',...
    'pos',[.65 .05 .325 .08],'callback',@showpop);
conpanel.uicontrols.overview = uicontrol('parent',PopLatFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .15 .25 .08],'callback',@LoadOverview);
conpanel.selectedidx = [];
conpanel.bw = [.5 .5 .5];
conpanel.rg = [1 0 0];

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(PopLatFig.popanel,'userdata',popanel)
set(gcf,'userdata',PopLatFig)


% Unpack Population data
UnpackPopData()
showpop()

% Display table
DispTable()


end

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end


function highlightcells(a,~)

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
popanel = get(PopLatFig.popanel,'userdata');
dispanel = get(PopLatFig.dispanel,'userdata');

% Which axes were selected?
if strcmp(a.Title.String,'Population Regressions (1D)')
    ndL = conpanel.oneDL;
else 
    ndL = conpanel.twoDL;
end

% Which latency measurement
if conpanel.uicontrols.latanalselector.Value == 1
    colsigL = conpanel.col.FSP.sigL;
    lumsigL = conpanel.lum.FSP.sigL;
else
    colsigL = conpanel.col.GTF.sigL;
    lumsigL = conpanel.lum.GTF.sigL;
end

% Index into subunits
lL = ndL & lumsigL;
cL = ndL & colsigL;


% Set colors on table
cols = repmat([1 1 1],[numel(lL) 1]);
cols(lL,:) = repmat(conpanel.bw,[sum(lL) 1]);
cols(cL,:) = repmat(conpanel.rg,[sum(cL) 1]);
cols(lL & cL,:) = repmat([1 .5 0],[sum(lL&cL) 1]);
%set(conpanel.table,'CellSelection;Callback',[])
set(conpanel.table,'BackgroundColor',cols)
%set(conpanel.table,'CellSelectionCallback',@cellselect

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(PopLatFig.popanel,'userdata',popanel)
set(gcf,'userdata',PopLatFig)


end

function cellselect(~,b)
global GLMSPopData

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
dispanel = get(PopLatFig.dispanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if ~isempty(b.Indices)
    conpanel.selectedidx = b.Indices(1);
end

% Plot the individual subunits
datatypes = GLMSPopData(1,:);
spikestats = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMPSpikeStats')};
colL = abs(spikestats.colprojcc) > .02;
lumL = abs(spikestats.lumprojcc) > .02;

% slopes and intercepts
colFSPint = conpanel.col.FSP.intercept(conpanel.selectedidx);
colFSPslo = conpanel.col.FSP.slope(conpanel.selectedidx);
colGTFint = conpanel.col.GTF.intercept(conpanel.selectedidx);
colGTFslo = conpanel.col.GTF.slope(conpanel.selectedidx);
lumFSPint = conpanel.lum.FSP.intercept(conpanel.selectedidx);
lumFSPslo = conpanel.lum.FSP.slope(conpanel.selectedidx);
lumGTFint = conpanel.lum.GTF.intercept(conpanel.selectedidx);
lumGTFslo = conpanel.lum.GTF.slope(conpanel.selectedidx);

% latencies and responses
respFSP = spikestats.firstsp.response;
latFSP = spikestats.firstsp.latency;
respGTF = spikestats.gauss.response;
latGTF = spikestats.gauss.latency;

xax = linspace(0,max([respFSP; respGTF]),50);

yaxcolFSP = xax*colFSPslo + colFSPint;
yaxcolGTF = xax*colGTFslo + colGTFint;
yaxlumFSP = xax*lumFSPslo + lumFSPint;
yaxlumGTF = xax*lumGTFslo + lumGTFint;

% Plot FSP
axes(dispanel.axes.FSP); cla; hold on; grid on;
plot(respFSP(lumL),latFSP(lumL),'o','color',conpanel.bw)
plot(respFSP(colL),latFSP(colL),'o','color',conpanel.rg)
plot(xax,yaxlumFSP,'--','color',conpanel.bw);
plot(xax,yaxcolFSP,'--','color',conpanel.rg);
xlabel('Response')
ylabel('Latency (s)')

% Plot GTF
axes(dispanel.axes.GTF); cla; hold on; grid on;
plot(respGTF(lumL),latGTF(lumL),'o','color',conpanel.bw)
plot(respGTF(colL),latGTF(colL),'o','color',conpanel.rg)
plot(xax,yaxlumGTF,'--','color',conpanel.bw);
plot(xax,yaxcolGTF,'--','color',conpanel.rg);
xlabel('Response')
ylabel('Latency (s)')

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(gcf,'userdata',PopLatFig)

showpop()


end

function callanalysis(~,~)
global GLMSPopData GLMP

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');

% If no cell selected, return prompt. If selected, call analysis.
if isempty(conpanel.selectedidx)
    disp('Must select a cell first')
else
    datatype = GLMSPopData(1,:);
    GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatype,'GLMP')};
    GLMSGUI_GLMPSpikeStats
end

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',PopLatFig)

end


function UnpackPopData()
global GLMSPopData

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
datatypes = GLMSPopData(1,:);

% unpack latency data
spikestats = GLMSPopData(2:end,strcmp(datatypes,'GLMPSpikeStats'));
h0 = [spikestats{:}];
h1 = [h0.regress];
firstsp = [h1.firstsp]';
gausssp = [h1.gauss]';
lumvalsFSP = [firstsp.lum];
colvalsFSP = [firstsp.col];
lumvalsGTF = [gausssp.lum];
colvalsGTF = [gausssp.col];
lumFSPcoefs = cat(2,lumvalsFSP.coefs)';% 1st col is offset, 2nd is slope
lumGTFcoefs = cat(2,lumvalsGTF.coefs)';
colFSPcoefs = cat(2,colvalsFSP.coefs)';
colGTFcoefs = cat(2,colvalsGTF.coefs)';
conpanel.lum.FSP.slope = lumFSPcoefs(:,2);
conpanel.lum.FSP.intercept = lumFSPcoefs(:,1);
conpanel.lum.FSP.pval = [lumvalsFSP.pval]';
conpanel.lum.GTF.slope = lumGTFcoefs(:,2);
conpanel.lum.GTF.intercept = lumGTFcoefs(:,1);
conpanel.lum.GTF.pval = [lumvalsGTF.pval]';
conpanel.col.FSP.slope = colFSPcoefs(:,2);
conpanel.col.FSP.intercept = colFSPcoefs(:,1);
conpanel.col.FSP.pval = [colvalsFSP.pval]';
conpanel.col.GTF.slope = colGTFcoefs(:,2);
conpanel.col.GTF.intercept = colGTFcoefs(:,1);
conpanel.col.GTF.pval = [colvalsGTF.pval]';

% nd params
nd = GLMSPopData(2:end,strcmp(datatypes,'nD'));
conpanel.oneDL = strcmp(nd,'1D');
conpanel.twoDL = strcmp(nd,'2D');
conpanel.nd = zeros(size(conpanel.oneDL));
conpanel.nd(conpanel.oneDL) = 1;
conpanel.nd(conpanel.twoDL) = 2;

% tuning params
tuning = [GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]';
degpm = 22;
conpanel.col.L = (tuning > -45-degpm & tuning < -45+degpm)...
    | (tuning > 135-degpm & tuning < 135+degpm);
conpanel.lum.L = (tuning > -135-degpm & tuning < -135+degpm)...
    | (tuning > 45-degpm & tuning < 45+degpm);

% Make a logical index for significance (pval)
sigval = 2;
conpanel.lum.FSP.sigL = conpanel.lum.FSP.pval < sigval;
conpanel.lum.GTF.sigL = conpanel.lum.GTF.pval < sigval;
conpanel.col.FSP.sigL = conpanel.col.FSP.pval < sigval;
conpanel.col.GTF.sigL = conpanel.col.GTF.pval < sigval;

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',PopLatFig)

end


function showpop(~,~)

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
dispanel = get(PopLatFig.dispanel,'userdata');
popanel = get(PopLatFig.popanel,'userdata');

% Which measure of latency?
if conpanel.uicontrols.latanalselector.Value == 1
    lumslos = conpanel.lum.FSP.slope;
    lumints = conpanel.lum.FSP.intercept;
    lumsigL = conpanel.lum.FSP.sigL;
    colslos = conpanel.col.FSP.slope;
    colints = conpanel.col.FSP.intercept;
    colsigL = conpanel.col.FSP.sigL;
else
    lumslos = conpanel.lum.GTF.slope;
    lumints = conpanel.lum.GTF.intercept;
    lumsigL = conpanel.lum.GTF.sigL;
    colslos = conpanel.col.GTF.slope;
    colints = conpanel.col.GTF.intercept;
    colsigL = conpanel.col.GTF.sigL;
end

% Organize by direction and nD (1D)
lL = conpanel.oneDL & lumsigL;
cL = conpanel.oneDL & colsigL;

% plot figure
axes(popanel.axes.oneD); cla; hold on;
plot(lumints(lL),lumslos(lL),'o','color',conpanel.bw);
plot(colints(cL),colslos(cL),'o','color',conpanel.rg);
xlabel('Intercept of latency on response')
ylabel('Slope of latency on response')
title('Population Regressions (1D)')

% Organize by direction and nD (2D)
lL = conpanel.twoDL & lumsigL;
cL = conpanel.twoDL & colsigL;

% plot figure
axes(popanel.axes.twoD); cla; hold on;
plot(lumints(lL),lumslos(lL),'o','color',conpanel.bw);
plot(colints(cL),colslos(cL),'o','color',conpanel.rg);
xlabel('Intercept of latency on response')
ylabel('Slope of latency on response')
title('Population Regressions (2D)')

% Set colors on table
%cols = repmat([1 1 1],[numel(lL) 1]);
%set(conpanel.table,'BackgroundColor',cols)

if ~isempty(conpanel.selectedidx)
    if conpanel.nd(conpanel.selectedidx) == 1
        axes(popanel.axes.oneD); hold on;
    else
        axes(popanel.axes.twoD); hold on;
    end
    plot(lumints(conpanel.selectedidx),lumslos(conpanel.selectedidx),'*k')
    plot(colints(conpanel.selectedidx),colslos(conpanel.selectedidx),'*k')
end


% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(gcf,'userdata',PopLatFig)

end


function DispTable()
global GLMSPopData

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
dispanel = get(PopLatFig.dispanel,'userdata');

% Unpack data
datatypes = GLMSPopData(1,:);
data = cell(size(GLMSPopData,1)-1,3);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Tuning'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'nD'));

% Repackage data and name columns
colname = {'Datafile' 'Tuning' 'nD'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table
set(conpanel.table,'data',data,'columnname',colname,...
    'columnformat',colformat);

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(gcf,'userdata',PopLatFig)

end

function Reanal(~,~,~)
global GLMSPopData GLMP

% Load figure and pop variables
PopLatFig = get(gcf,'userdata');
conpanel = get(PopLatFig.conpanel,'userdata');
dispanel = get(PopLatFig.dispanel,'userdata');

% Set up variables
datatypes = GLMSPopData(1,:);
subs = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
GLMPs = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Run through each dataset
for n = 1:numel(GLMPs)
    GLMP = GLMPs{n};
    sub = subs{n};
    GLMSGUI_GLMPSpikeStats(sub)
    spikefig = get(7548,'userdata');
    scatfig = spikefig.scatfig;
    GLMSPopData{n+1,strcmp(datatypes,'GLMPSpikeStats')} = scatfig;

    % Save population data at each iteration
    save([library 'GLMSPopData'],'GLMSPopData');
    
end

disp('Population Latency Analysis Completed')

% Save user data
set(PopLatFig.conpanel,'userdata',conpanel)
set(PopLatFig.dispanel,'userdata',dispanel)
set(gcf,'userdata',PopLatFig)

end
