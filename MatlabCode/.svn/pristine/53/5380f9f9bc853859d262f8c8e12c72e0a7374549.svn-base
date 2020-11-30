function GLMSPopGUI_Regress(~,~)
% For plotting within the GLMSPop Regression Analysis
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

% Set up figure and panels
figure(6447); clf; 
set(gcf,'NumberTitle','off','Name','Population Regression','units','normalized',...
    'pos',[.1 .1 .7 .7])
regressfig.table = uitable('units','normalized','pos',[.01 .2 .24 .79],...
    'BackgroundColor',[1 1 1]);%,'cellselectioncallback',@cellselect);
regressfig.xvarpanel = uibuttongroup('units','normalized','pos',[.26 .6 .13 .39],...
    'title','X Variable','SelectionChangedFcn',@RegressPlot);
regressfig.yvarpanel = uibuttongroup('units','normalized','pos',[.26 .2 .13 .39],...
    'title','Y Variable','SelectionChangedFcn',@RegressPlot);
regressfig.conpanel = uipanel('units','normalized','pos',[.01 .01 .38 .18]);
regressfig.plot = axes('parent',gcf,'units','normalized','pos',[.5 .1 .4 .8]); hold on; grid on; box on;

% Populate table
datatypes = GLMSPopData(1,:);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'nD'));
data(:,4) = num2cell(([GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]'));
colname = {'Datafile' 'Sub' 'nD' 'Tuning'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
set(regressfig.table,'data',data,'columnname',colname,...
    'columnformat',colformat,'CellSelectionCallback',@cellselect);

% Set up variable panels
strs = {'xvarpanel' 'yvarpanel'};
for n = 1:2
    vars.(strs{n}).nstim = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .9 .8 .05],'style','radiobutton','String','# Stimuli');
    vars.(strs{n}).prefdir = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .8 .8 .05],'style','radiobutton','String','Pref Dir');
    vars.(strs{n}).prefdirvar = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .7 .8 .05],'style','radiobutton','String','Hessian PD Std');
    vars.(strs{n}).prefdirvar = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .6 .8 .05],'style','radiobutton','String','BF PD Std');
    vars.(strs{n}).blfr = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .5 .8 .05],'style','radiobutton','String','BL FR');
    vars.(strs{n}).kappa = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .4 .8 .05],'style','radiobutton','String','Kappa');
    vars.(strs{n}).topresp = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .3 .8 .05],'style','radiobutton','String','Top Mean Resp');
    vars.(strs{n}).date = uicontrol('parent',regressfig.(strs{n}),'units','normalized',...
        'pos',[.1 .2 .8 .05],'style','radiobutton','String','Date');
end

% Preallocate space for variables
conpanel.vars.nstim = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.subs = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.prefdir = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.prefdirstdHess = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.prefdirstdBF = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.nd = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.blfr = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.kappa = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.topmeanresp = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.date = (1:size(GLMSPopData,1)-1)';
conpanel.selectedidx = [];

% Pull out the variables to be plotted
datatypes = GLMSPopData(1,:);
datafiles = cat(1,GLMSPopData{2:end,strcmp(datatypes,'Datafile')});
GLMPs = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
conpanel.vars.subs = [GLMSPopData{2:end,strcmp(datatypes,'Subunit')}]';
conpanel.vars.prefdir = [GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]';
nds = cat(1,GLMSPopData{2:end,strcmp(datatypes,'nD')});
surfs = GLMSPopData(2:end,strcmp(datatypes,'Surface Parameters'));
tunstd = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';

strs = {'oneD' 'twoD'};
for n = 1:numel(GLMPs)
    conpanel.vars.nd(n) = str2double(nds(n,1));
    conpanel.vars.nstim(n) = numel(GLMPs{n}.subunit{conpanel.vars.subs(n)}.Lcc);
    conpanel.vars.prefdirstdHess(n) = sqrt(1./surfs{n}.(strs{conpanel.vars.nd(n)}).Hessian(end-1,end-1))/pi*180;
    conpanel.vars.prefdirstdBF(n) = tunstd(n).bfconfdeg.(strs{conpanel.vars.nd(n)});
    conpanel.vars.blfr(n) = surfs{n}.(strs{conpanel.vars.nd(n)}).parvals(end-2);
    conpanel.vars.kappa(n) = surfs{n}.(strs{conpanel.vars.nd(n)}).parvals(end);
    conpanel.vars.topmeanresp(n) = max(GLMPs{n}.subunit{conpanel.vars.subs(n)}.meannspikes);
end

% Set up controls for conpanel
conpanel.uicontrols.overview = uicontrol('parent',regressfig.conpanel,...
    'units','normalized','pos',[.05 .1 .425 .35],...
    'style','pushbutton','string','Load Overview','callback',@LoadOverview);

% Set up selector for which monkey and nd
strs = {'All Data' 'Nut Data' 'Maui Data'};
conpanel.uicontrols.whichMonk = uicontrol('parent',regressfig.conpanel,...
    'style','popup','units','normalized','pos',[.6 .2 .3 .2],...
    'string',strs,'value',1,'callback',@RegressPlot);
strs = {'1D & 2D' '1D' '2D'};
conpanel.uicontrols.whichnD = uicontrol(regressfig.conpanel,...
    'style','popup','units','normalized','pos',[.6 .6 .3 .2],...
    'string',strs,'value',1,'callback',@RegressPlot);

% Set up monkey index
conpanel.monkL.Nut = datafiles(:,1) == 'N';
conpanel.monkL.Maui = datafiles(:,1) == 'M';
conpanel.monkL.all = ones(size(datafiles,1),1);

% Set up nd idx
conpanel.ndL.oneD = conpanel.vars.nd == 1;
conpanel.ndL.twoD = conpanel.vars.nd == 2;
conpanel.ndL.all = conpanel.monkL.all;

% Save figure variables
set(regressfig.conpanel,'userdata',conpanel)
set(regressfig.xvarpanel,'userdata',vars.xvarpanel)
set(regressfig.yvarpanel,'userdata',vars.yvarpanel)
%set(regressfig.plot,'userdata',plotpan)
set(gcf,'userdata',regressfig)

% Call plotting function
RegressPlot

end

function cellselect(~,b)

% Load figure and pop variables
regressfig = get(6447,'userdata');
conpanel = get(regressfig.conpanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if ~isempty(b.Indices)
    conpanel.selectedidx = b.Indices(1);
else
    return
end

% Save figure variables
set(regressfig.conpanel,'userdata',conpanel)
set(6447,'userdata',regressfig)

RegressPlot()

end



function RegressPlot(~,~)

%disp('regress plot called')

% Unpack figure variables
regressfig = get(gcf,'userdata');
conpanel = get(regressfig.conpanel,'userdata');
%xvarpan = get(regressfig.xvarpanel,'userdata');
%yvarpan = get(regressfig.yvarpanel,'userdata');
plotpan = get(regressfig.plot,'userdata');

% Finding selected data
strs = {'xvarpanel' 'yvarpanel'};
%vars = cell(1,2);
labels = cell(1,2);
for n = 1:2
    str = regressfig.(strs{n}).SelectedObject.String;
    if strcmp(str,'# Stimuli')
        plotpan.vars{n} = conpanel.vars.nstim;
        labels{n} = str;
    elseif strcmp(str,'Pref Dir')
        plotpan.vars{n} = conpanel.vars.prefdir;
        labels{n} = 'Preferred Direction (deg)';
    elseif strcmp(str,'Hessian PD Std')
        plotpan.vars{n} = conpanel.vars.prefdirstdHess;
        labels{n} = 'Hessian PD Std (deg)';
    elseif strcmp(str,'BF PD Std')
        plotpan.vars{n} = conpanel.vars.prefdirstdBF;
        labels{n} = 'BF PD Std (deg)';  
    elseif strcmp(str,'BL FR')
        plotpan.vars{n} = conpanel.vars.blfr;
        labels{n} = 'Baseline (spike count)';
    elseif strcmp(str,'Kappa')
        plotpan.vars{n} = conpanel.vars.kappa;
        labels{n} = 'Kappa Value';
    elseif strcmp(str,'Top Mean Resp')
        plotpan.vars{n} = conpanel.vars.topmeanresp;
        labels{n} = 'Maximum Mean Response (spike count)';
    elseif strcmp(str,'Date')
        plotpan.vars{n} = conpanel.vars.date;
        labels{n} = 'Date (by rank order)';
    end
end

% Monkey idx
if strcmp(conpanel.uicontrols.whichMonk.String{conpanel.uicontrols.whichMonk.Value},'All Data')
    monkL = conpanel.monkL.all;
elseif strcmp(conpanel.uicontrols.whichMonk.String{conpanel.uicontrols.whichMonk.Value},'Nut Data')
    monkL = conpanel.monkL.Nut;
elseif strcmp(conpanel.uicontrols.whichMonk.String{conpanel.uicontrols.whichMonk.Value},'Maui Data')
    monkL = conpanel.monkL.Maui;
end

% nD idx
if strcmp(conpanel.uicontrols.whichnD.String{conpanel.uicontrols.whichnD.Value},'1D & 2D')
    ndL = conpanel.ndL.all;
elseif strcmp(conpanel.uicontrols.whichnD.String{conpanel.uicontrols.whichnD.Value},'1D')
    ndL = conpanel.ndL.oneD;
elseif strcmp(conpanel.uicontrols.whichnD.String{conpanel.uicontrols.whichnD.Value},'2D')
    ndL = conpanel.ndL.twoD;
end

% Plot data
plotpan.L = monkL & ndL;
axes(regressfig.plot); cla; hold on; grid on;
plot(plotpan.vars{1}(plotpan.L),plotpan.vars{2}(plotpan.L),'ko','ButtonDownFcn',@SelectPts)
xlabel(labels{1})
ylabel(labels{2})

if ~isempty(conpanel.selectedidx)
    if any(conpanel.selectedidx == find(plotpan.L))
        
        % Highlight data in plot
        plot(plotpan.vars{1}(conpanel.selectedidx),plotpan.vars{2}(conpanel.selectedidx),'*','color',[1 .5 0])
        
        % Highlight cells in table
        cols = repmat([1 1 1],[size(regressfig.table.Data,1) 1]);
        cols(conpanel.selectedidx,:) = [1 .5 0];
        %set(regressfig.table,'CellSelectionCallback',[])
        %keyboard
        set(regressfig.table,'BackgroundColor',cols)
        %set(regressfig.table,'CellSelectionCallback',@cellselect)
        
    else
        disp('NOTE: Selected Cell does not meet plotting criteria...')
        
        % Plot data even though it doesn't meet criteria
        plot(plotpan.vars{1}(conpanel.selectedidx),plotpan.vars{2}(conpanel.selectedidx),'ro','color',[1 0 0])
 
        % Highlight cells in table
        cols = repmat([1 1 1],[size(regressfig.table.Data,1) 1]);
        cols(conpanel.selectedidx,:) = [1 0 0];
        set(regressfig.table,'CellSelectionCallback',[])
        set(regressfig.table,'BackgroundColor',cols)
    end
end
%set(regressfig.table,'CellSelectionCallback',@cellselect);

% Save figure variables
set(regressfig.conpanel,'userdata',conpanel)
%set(regressfig.xvarpanel,'userdata',vars.xvarpanel)
%set(regressfig.yvarpanel,'userdata',vars.yvarpanel)
set(regressfig.plot,'userdata',plotpan)
set(gcf,'userdata',regressfig)

end

function SelectPts(~,~)
%disp('selectpts called')

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Unpack figure variables
regressfig = get(gcf,'userdata');
conpanel = get(regressfig.conpanel,'userdata');
%xvarpan = get(regressfig.xvarpanel,'userdata');
%yvarpan = get(regressfig.yvarpanel,'userdata');
plotpan = get(regressfig.plot,'userdata');

% Which point was selected? Find out!
[~,c] = min(sqrt((plotpan.vars{1}(plotpan.L)-whichpt(1)).^2 + (plotpan.vars{2}(plotpan.L)-whichpt(2)).^2));
tmp = find(plotpan.L);
conpanel.selectedidx = tmp(c);

% Save figure variables
set(regressfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',regressfig)

RegressPlot();

end

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
regressfig = get(gcf,'userdata');
conpanel = get(regressfig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

