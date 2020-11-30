function GLMSPopGUI_SelectCells
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])


% Set up figure and panels
figure(6448); clf; 
set(gcf,'NumberTitle','off','Name','Population Cell Selection','units','normalized',...
    'pos',[.1 .1 .8 .8])
celselfig.table = uitable('units','normalized','pos',[.01 .51 .48 .48],...
    'BackgroundColor',[1 1 1]);
celselfig.conpanel = uipanel('units','normalized','pos',[.01 .01 .48 .48]);
celselfig.disp = axes('units','normalized','pos',[.55 .1 .4 .8]); hold on; grid on; box on;

%%% Control Panel %%%

% Std of Pref Dir
conpanel.PDstdpan = uipanel('parent',celselfig.conpanel,'units','normalized',...
    'pos',[.01 .76 .98 .2],'title','Std of Preferred Direction (deg)');
vars.PDstd.min = uicontrol('parent',conpanel.PDstdpan,'style','edit',...
    'units','normalized','Pos',[.05 .1 .1 .8],'String','0',...
    'callback',@PlotData);
vars.PDstd.max = uicontrol('parent',conpanel.PDstdpan,'style','edit',...
    'units','normalized','Pos',[.2 .1 .1 .8],'String','Inf',...
    'callback',@PlotData);
vars.PDstd.histax = axes('parent',conpanel.PDstdpan,'units','normalize',...
    'pos',[.4 .3 .55 .6],'box','on');

% Modulation
conpanel.modpan = uipanel('parent',celselfig.conpanel,'units','normalized',...
    'pos',[.01 .52 .98 .2],'title','Modulation from Baseline (sp/s)');
vars.mod.min = uicontrol('parent',conpanel.modpan,'style','edit',...
    'units','normalized','Pos',[.05 .1 .1 .8],'String','0',...
    'callback',@PlotData);
vars.mod.max = uicontrol('parent',conpanel.modpan,'style','edit',...
    'units','normalized','Pos',[.2 .1 .1 .8],'String','Inf',...
    'callback',@PlotData);
vars.mod.histax = axes('parent',conpanel.modpan,'units','normalized',...
    'pos',[.4 .3 .55 .6],'box','on');

% Baseline
conpanel.BLpan = uipanel('parent',celselfig.conpanel,'units','normalized',...
    'pos',[.01 .28 .98 .2],'title','Mean Basline Response (sp/s)');
vars.BL.min = uicontrol('style','edit','parent',conpanel.BLpan,...
    'units','normalized','Pos',[.05 .1 .1 .8],'String','0',...
    'callback',@PlotData);
vars.BL.max = uicontrol('style','edit','parent',conpanel.BLpan,...
    'units','normalized','Pos',[.2 .1 .1 .8],'String','Inf',...
    'callback',@PlotData);
vars.BL.histax = axes('parent',conpanel.BLpan,'units','normalized',...
    'pos',[.4 .3 .55 .6],'box','on');

% Kappa
conpanel.kappapan = uipanel('parent',celselfig.conpanel,'units','normalized',...
    'pos',[.01 .04 .98 .2],'title','Kappa');
vars.kappa.min = uicontrol('parent',conpanel.kappapan,'style','edit',...
    'units','normalized','Pos',[.05 .1 .1 .8],'String','0',...
    'callback',@PlotData);
vars.kappa.max = uicontrol('parent',conpanel.kappapan,'style','edit',...
    'units','normalized','Pos',[.2 .1 .1 .8],'String','Inf',...
    'callback',@PlotData);
vars.kappa.histax = axes('parent',conpanel.kappapan,'units','normalize',...
    'pos',[.4 .3 .55 .6],'box','on');

%Save variables
conpanel.idx = 0;
conpanel.vars = vars;
set(celselfig.conpanel,'UserData',conpanel);
set(gcf,'UserData',celselfig);

UnpackPopData()

PlotData()


end


function PlotData(~,~)

% Load user data
celselfig = get(gcf,'userdata');
conpanel = get(celselfig.conpanel,'UserData');

% Find which points are within the range
PDstdL = str2double(conpanel.vars.PDstd.min.String) <= conpanel.vars.PDstd.vals...
    & conpanel.vars.PDstd.vals <= str2double(conpanel.vars.PDstd.max.String);
modL = str2double(conpanel.vars.mod.min.String) <= conpanel.vars.mod.vals...
    & conpanel.vars.mod.vals <= str2double(conpanel.vars.mod.max.String);
BLL = str2double(conpanel.vars.BL.min.String) <= conpanel.vars.BL.vals...
    & conpanel.vars.BL.vals <= str2double(conpanel.vars.BL.max.String);
kappaL = str2double(conpanel.vars.kappa.min.String) <= conpanel.vars.kappa.vals...
    & conpanel.vars.kappa.vals <= str2double(conpanel.vars.kappa.max.String);
goodL = PDstdL & modL & BLL & kappaL;

% Plot good subunits and bad subunits
axes(celselfig.disp); cla; hold on; grid on; box on;
h1 = plot3(conpanel.vars.PDstd.vals(goodL),conpanel.vars.mod.vals(goodL),...
    conpanel.vars.BL.vals(goodL),'ko','MarkerFaceColor','g');
h2 = plot3(conpanel.vars.PDstd.vals(~goodL),conpanel.vars.mod.vals(~goodL),...
    conpanel.vars.BL.vals(~goodL),'ko','MarkerFaceColor','r');
if conpanel.idx ~= 0
    plot3(conpanel.vars.PDstd.vals(conpanel.idx),conpanel.vars.mod.vals(conpanel.idx),...
        conpanel.vars.BL.vals(conpanel.idx),'Marker','pentagram',...
        'MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b');
end
xlabel('PDstd');
ylabel('mod')
zlabel('BL')
set(celselfig.disp,'ButtonDownFcn',@cellselect);
set(h1,'ButtonDownFcn',@cellselect);
set(h2,'ButtonDownFcn',@cellselect);

%%% Plot Histograms %%%
% PDstd
bins = linspace(0,max(conpanel.vars.PDstd.vals),50);
axes(conpanel.vars.PDstd.histax); cla; hold on;
[a,b] = hist(conpanel.vars.PDstd.vals,bins);
bar(b,a)
xlim([bins(1) bins(end)])
ylim([0 max(a)])
plot(repmat(str2double(conpanel.vars.PDstd.min.String),[1 2]),[min(a) max(a)],'r--')
plot(repmat(str2double(conpanel.vars.PDstd.max.String),[1 2]),[min(a) max(a)],'r--')

% Mod
bins = linspace(0,max(conpanel.vars.mod.vals),50);
axes(conpanel.vars.mod.histax); cla; hold on;
[a,b] = hist(conpanel.vars.mod.vals,bins);
bar(b,a)
xlim([bins(1) bins(end)])
ylim([0 max(a)])
plot(repmat(str2double(conpanel.vars.mod.min.String),[1 2]),[min(a) max(a)],'r--')
plot(repmat(str2double(conpanel.vars.mod.max.String),[1 2]),[min(a) max(a)],'r--')

% BL
bins = linspace(0,max(conpanel.vars.BL.vals),50);
axes(conpanel.vars.BL.histax); cla; hold on;
[a,b] = hist(conpanel.vars.BL.vals,bins);
bar(b,a)
xlim([bins(1) bins(end)])
ylim([0 max(a)])
plot(repmat(str2double(conpanel.vars.BL.min.String),[1 2]),[min(a) max(a)],'r--')
plot(repmat(str2double(conpanel.vars.BL.max.String),[1 2]),[min(a) max(a)],'r--')

% Kappa
bins = linspace(0,max(conpanel.vars.kappa.vals),50);
axes(conpanel.vars.kappa.histax); cla; hold on;
[a,b] = hist(conpanel.vars.kappa.vals,bins);
bar(b,a)
xlim([bins(1) bins(end)])
ylim([0 max(a)])
plot(repmat(str2double(conpanel.vars.kappa.min.String),[1 2]),[min(a) max(a)],'r--')
plot(repmat(str2double(conpanel.vars.kappa.max.String),[1 2]),[min(a) max(a)],'r--')

% Set table colors
PopulateTable(goodL)

% Save user variables
conpanel.goodL = goodL;
set(celselfig.conpanel,'UserData',conpanel);
set(gcf,'UserData',celselfig);

end

function PopulateTable(varargin)
global GLMSPopData

% Load user data
celselfig = get(gcf,'userdata');
conpanel = get(celselfig.conpanel,'UserData');

% Populate table
datatypes = GLMSPopData(1,:);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'nD'));
data(:,4) = num2cell(([GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]'));
data(:,5) = num2cell(conpanel.vars.PDstd.vals);
data(:,6) = num2cell(conpanel.vars.mod.vals);
data(:,7) = num2cell(conpanel.vars.BL.vals);
data(:,8) = num2cell(conpanel.vars.kappa.vals);
colname = {'Datafile' 'Sub' 'nD' 'Tuning' 'Tun Std' 'Mod' 'BL' 'Kappa'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
cols = ones(size(data,1),3);
if ~isempty(varargin)
    cols(~varargin{1},:) = repmat([1 0 0],sum(~varargin{1}),1);
end
if conpanel.idx ~= 0
    cols(conpanel.idx,:) = [.5 .5 1];
end
set(celselfig.table,'data',data,'columnname',colname,'BackgroundColor',cols,...
    'columnformat',colformat,'busyaction','cancel','Interruptible','off',...
    'CellSelectionCallback',@cellselect);

end

function cellselect(a,b)

% Load user data
celselfig = get(gcf,'userdata');
conpanel = get(celselfig.conpanel,'UserData');
if strcmp(get(a,'type'),'uitable')
    pause(.5) % This is totally stupid but necessary to not re-select the previous cell
    if isempty(b.Indices) || b.Indices(1) == conpanel.idx
        return
    else
       conpanel.idx = b.Indices(1); 
    end
elseif strcmp(get(a,'type'),'line') || strcmp(get(a,'type'),'axes')
    % Find dataset with minimum Euclidean distance from selected point
    pdstddist = conpanel.vars.PDstd.vals - b.IntersectionPoint(1);
    moddist = conpanel.vars.mod.vals - b.IntersectionPoint(2);
    bldist = conpanel.vars.BL.vals - b.IntersectionPoint(3);
    [~,conpanel.idx] = min(sqrt(pdstddist.^2 + moddist.^2 + bldist.^2));
end

% Highlight cell
% if ~isempty(varargin)
%     cols(~varargin{1},:) = repmat([1 0 0],sum(~varargin{1}),1);
% end
% if conpanel.idx ~= 0
%     cols(conpanel.idx,:) = [.5 .5 1];
% end

%Save variables
set(celselfig.table,'enable','off')
set(celselfig.conpanel,'UserData',conpanel);
set(gcf,'UserData',celselfig);

PlotData();
set(celselfig.table,'enable','on')

end

function UnpackPopData()
global GLMSPopData

% Load user data
celselfig = get(gcf,'userdata');
conpanel = get(celselfig.conpanel,'UserData');

% Load data
datatypes = GLMSPopData(1,:);

% Baseline and Modulation
conpanel.vars.BL.vals = nan(size(GLMSPopData,1)-1,1);
conpanel.vars.mod.vals = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    conpanel.vars.BL.vals(n-1) = mean(GLMP.subunit{sub}.blfr);
    conpanel.vars.mod.vals(n-1) = max(GLMP.subunit{sub}.meannspikes)...
        ./ mean(GLMP.subunit{sub}.stimDur) - mean(GLMP.subunit{sub}.blfr);
end

% Pull out population parameters
nds = GLMSPopData(2:end,strcmp(datatypes,'nD'));
onedL = strcmp(nds,'1D');
poparams = [GLMSPopData{2:end,strcmp(datatypes,'Surface Parameters')}];
params = nan(numel(nds),numel(poparams(1).oneD.parvals));
onedpars = [poparams.oneD];
twodpars = [poparams.twoD];
temp1 = cat(1,onedpars.parvals);
temp2 = cat(1,twodpars.parvals);
params(onedL,:) = temp1(onedL,:);
params(~onedL,:) = temp2(~onedL,:);

% Kappa
conpanel.vars.kappa.vals = params(:,end);

% Tuning Std
conpanel.vars.PDstd.vals = nan(size(conpanel.vars.kappa.vals));
confint = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';
t = [confint.bfconfdeg]';
onedconf = [t.oneD]';
twodconf = [t.twoD]';
conpanel.vars.PDstd.vals(onedL) = onedconf(onedL);
conpanel.vars.PDstd.vals(~onedL) = twodconf(~onedL);

%Save variables
set(celselfig.conpanel,'UserData',conpanel);
set(gcf,'UserData',celselfig);

end



