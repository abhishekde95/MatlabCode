function GLMSPopGUI_RFLocations(~,~)
global GLMSPopData

figure(63); clf;
set (gcf,'pos',[200 200 800 500],'numbertitle','off','name','Receptive Field Locations Analysis');
RFFig.conpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.01 .01 .34 .38]);
RFFig.tuningpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.01 .4 .34 .59]);
RFFig.RFpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.36 .01 .68 .98]);

% Look at 1D or 2D cells
strs = {'1D Cells' '2D Cells' 'All Cells'};
conpanel.nDmenu = uicontrol(RFFig.conpanel,'style','popupmenu',...
    'units','normalized','pos',[.05 .8 .4 .15],'string',strs,'Callback',@resetPlots);

% Look at all directions or selected directions
strs = {'All Directions' 'Compare Directions'};
conpanel.tunedirmenu = uicontrol(RFFig.conpanel,'style','popupmenu',...
    'units','normalized','pos',[.55 .8 .4 .15],'string',strs,'Callback',@resetPlots);

% Buttons for comparing different tune directions
strs = {'Bin #1' 'Bin #2' 'Bin #3' 'Bin #4'};
conpanel.binselectormenu = uicontrol('parent',RFFig.conpanel,'style','popupmenu',...
    'units','normalized','pos',[.05 .6 .4 .15],...
    'string',strs,'enable','off');

% Monkey menu
strs = {'Nut Data' 'Maui Data' 'All Data'};
conpanel.monkeyselmenu = uicontrol('parent',RFFig.conpanel,'style','popupmenu',...
    'units','normalized','pos',[.55 .6 .4 .15],'string',strs,'Callback',@resetPlots);

conpanel.reanalall = uicontrol('parent',RFFig.conpanel,'style','pushbutton',...
    'units','normalized','pos',[.1 .1 .8 .2],'string','Reanalyze All','callback',@ReanalAll);

% Set up RF Axes
rfpanel.axes = axes('parent',RFFig.RFpanel,'units','normalized','pos',[.1 .1 .8 .8]);

% Set up Tuning axes
tpanel.axes = axes('parent',RFFig.tuningpanel,'unit','normalized',...
    'pos',[.1 .1 .8 .8]);

% Pull out some population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
datanames = GLMSPopData(1,:);
nds = GLMSPopData(2:end,strcmp(datanames,'nD'));
tpanel.oneDL = strcmp(nds,'1D');
tpanel.twoDL = strcmp(nds,'2D');
tpanel.allDL = true(size(nds));
datafiles = GLMSPopData(2:end,strcmp(datanames,'Datafile'));
datafiles = cat(1,datafiles{:});
tpanel.NutL = strcmp(cellstr(datafiles(:,1)),'N');
tpanel.MauiL = strcmp(cellstr(datafiles(:,1)),'M');
tpanel.allmonkeyL = true(size(tpanel.NutL));
tpanel.tuning = cat(1,GLMSPopData{2:end,strcmp(datanames,'Tuning')})./180*pi;
rfs = cat(1,GLMSPopData{2:end,strcmp(datanames,'RF Params')});
tpanel.rfthetas = [rfs.rftheta]';
tpanel.rfrhos = [rfs.rfrho]';

% Save User Data
set(RFFig.tuningpanel,'userdata',tpanel);
set(RFFig.RFpanel,'userdata',rfpanel);
set(RFFig.conpanel,'userdata',conpanel);
set(gcf,'userdata',RFFig);

PlotAllDirs()

end

function resetPlots(~,~)

% Load figure variables
RFFig = get(gcf,'userdata');
conpanel = get(RFFig.conpanel,'userdata');
tpanel = get(RFFig.tuningpanel,'userdata');
rfpanel = get(RFFig.RFpanel,'userdata');

if conpanel.tunedirmenu.Value == 1
    
    % Plot RFs
    PlotAllDirs();
    
    % Turn bin buttons off
    set(conpanel.binselectormenu,'enable','off')

elseif conpanel.tunedirmenu.Value == 2
    
    % Which nd
    if conpanel.nDmenu.Value == 1
        ndL = tpanel.oneDL;
    elseif conpanel.nDmenu.Value == 2
        ndL = tpanel.twoDL;
    elseif conpanel.nDmenu.Value == 3
        ndL = tpanel.allDL;
    end
    
    % which monkey
    if conpanel.monkeyselmenu.Value == 1
        monL = tpanel.NutL;
    elseif conpanel.monkeyselmenu.Value == 2
        monL = tpanel.MauiL;
    elseif conpanel.monkeyselmenu.Value == 3
        monL = tpanel.allmonkeyL;
    end
    
    % Plot tuning histogram
    axes(tpanel.axes); cla;
    tpanel.hist = rose(tpanel.tuning(ndL & monL));
    set(tpanel.hist,'ButtonDownFcn',@PlotCompareDirs)
    set(tpanel.axes,'ButtonDownFcn',@PlotCompareDirs)

    % Clear rf panel
    axes(rfpanel.axes); cla;

    % Turn on bin buttons
    set(conpanel.binselectormenu,'enable','on')
    for n = 1:4
        rfpanel.bin{n}.rfthetas = [];
        rfpanel.bin{n}.rfrhos = [];
        tpanel.bin{n}.angs = [];
    end
    set(tpanel.hist,'color','k','ButtonDownFcn',@PlotCompareDirs)
else
    keyboard
end

% Save figure variables
set(RFFig.tuningpanel,'userdata',tpanel)
set(RFFig.conpanel,'userdata',conpanel)
set(RFFig.RFpanel,'userdata',rfpanel)
set(gcf,'userdata',RFFig);

end

function PlotCompareDirs(~,~)

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);
[theta,~] = cart2pol(whichpt(1),whichpt(2));
theta = theta+pi;

% Load Variables
RFFig = get(gcf,'UserData');
tpanel = get(RFFig.tuningpanel,'UserData');
conpanel = get(RFFig.conpanel,'UserData');
rfpanel = get(RFFig.RFpanel,'userdata');

% Which nd
if conpanel.nDmenu.Value == 1
    ndL = tpanel.oneDL;
elseif conpanel.nDmenu.Value == 2
    ndL = tpanel.twoDL;
elseif conpanel.nDmenu.Value == 3
    ndL = tpanel.allDL;
end

% which monkey
if conpanel.monkeyselmenu.Value == 1
    monL = tpanel.NutL;
elseif conpanel.monkeyselmenu.Value == 2
    monL = tpanel.MauiL;
elseif conpanel.monkeyselmenu.Value == 3
    monL = tpanel.allmonkeyL;
end

% parameters
angs = tpanel.tuning(ndL & monL);
rfthetas = tpanel.rfthetas(ndL & monL);
rfrhos = tpanel.rfrhos(ndL & monL);

% Set up bin colors
cols = ['r';'g';'b';'m'];
whichbin = conpanel.binselectormenu.Value;

% Find the edges of the selected TuneDir
[tout,~] = rose(angs);
edges = unique(tout);
lowerTheta = edges(find(sign(theta-edges)==1,1,'last'))-pi;
upperTheta = edges(find(sign(theta-edges)==-1,1,'first'))-pi;

% Find and plot selected angles
angsIdx = find(angs >= lowerTheta & angs <= upperTheta);
rfpanel.bin{whichbin}.rfthetas = rfthetas(angsIdx);
rfpanel.bin{whichbin}.rfrhos = rfrhos(angsIdx);
tpanel.bin{whichbin}.angs = angs(angsIdx);

% Pull out greatest eccentricity (to keep graph consistent dimensions)
maxtheta = ceil(max(tpanel.rfrhos));

% Plot RF Locations
axes(rfpanel.axes); cla; hold on;
for n = 1:4
    polar(rfpanel.bin{n}.rfthetas,rfpanel.bin{n}.rfrhos,[cols(n) '*'])
end
set(rfpanel.axes,'xlim',[-maxtheta maxtheta],'ylim',[-maxtheta maxtheta])

% Plot Tuning Histogram
axes(tpanel.axes); cla; hold on;
for n = 1:4
    if ~isempty(tpanel.bin{n}.angs)
        rose2(tpanel.bin{n}.angs,[],cols(n))
    end
end
tpanel.hist = rose(angs); hold off;
set(tpanel.hist,'color','k','ButtonDownFcn',@PlotCompareDirs)

% Save variables
set(RFFig.RFpanel,'userdata',rfpanel);
set(RFFig.conpanel,'userdata',conpanel);
set(RFFig.tuningpanel,'userdata',tpanel);
set(gcf,'userdata',RFFig);

end

function PlotAllDirs(~,~)

% Load Variables
RFFig = get(gcf,'UserData');
tpanel = get(RFFig.tuningpanel,'UserData');
conpanel = get(RFFig.conpanel,'UserData');
rfpanel = get(RFFig.RFpanel,'userdata');

% Which nD
if conpanel.nDmenu.Value == 1
    ndL = tpanel.oneDL;
elseif conpanel.nDmenu.Value == 2
    ndL = tpanel.twoDL;
elseif conpanel.nDmenu.Value == 3
    ndL = tpanel.allDL;
end

% Which monkney
if conpanel.monkeyselmenu.Value == 1
    monL = tpanel.NutL;
elseif conpanel.monkeyselmenu.Value == 2
    monL = tpanel.MauiL;
elseif conpanel.monkeyselmenu.Value == 3
    monL = tpanel.allmonkeyL;
end

angs = tpanel.tuning(ndL & monL);
rfthetas = tpanel.rfthetas(ndL & monL);
rfrhos = tpanel.rfrhos(ndL & monL);

% Plot tuning histogram
axes(tpanel.axes);
tpanel.hist = rose2(angs,[],'b');

% Plot RF locations
%maxrho = max(tpanel.rfrhos);

% Plot RF locations
axes(rfpanel.axes); cla;
rfpanel.allbins.plot = polar(rfthetas,rfrhos,'k*');
axis equal
%set(rfpanel.axes,'xlim',[-maxrho maxrho])

% Save variables
set(RFFig.RFpanel,'userdata',rfpanel);
set(RFFig.conpanel,'userdata',conpanel);
set(RFFig.tuningpanel,'userdata',tpanel);
set(gcf,'userdata',RFFig);

end


function ReanalAll(~,~,~)
global GLMSPopData

% Load Variables
RFFig = get(gcf,'UserData');
tpanel = get(RFFig.tuningpanel,'UserData');
conpanel = get(RFFig.conpanel,'UserData');
rfpanel = get(RFFig.RFpanel,'userdata');

if ismac
    library1 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
    library2 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
else
    library1 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
    library2 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
datatypes = GLMSPopData(1,:);

disp('Re-calculating all RF parameters...');
for n = 2:(size(GLMSPopData,1))
    
    datafile = GLMSPopData{n,strcmp(datatypes,'Datafile')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    rawdata = nex2stro([library1 char(datafile) '.nex']);
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    DN = GLMSPopData{n,strcmp(datatypes,'DN')};
    GLMP.rf_x = rawdata.sum.exptParams.rf_x/10;
    GLMP.rf_y = rawdata.sum.exptParams.rf_y/10;
    [theta,rho] = cart2pol(GLMP.rf_x,GLMP.rf_y);
%     if rho > 7
%         rho
%         n
%         keyboard
%     end
    rf.rfx = GLMP.rf_x;
    rf.rfy = GLMP.rf_y;
    rf.rftheta = theta;
    rf.rfrho = rho;
    rf.dvaperstix = DN.DVAPerStix(1);
    rf.nstix = numel(GLMP.subunit{sub}.gridX{1});
    rf.sqdva = DN.DVAPerStix(1)*rf.nstix;
    
    
    % Save recalculated rf params
    GLMSPopData{n,strcmp(datatypes,'RF Params')} = rf;
    GLMSPopData{n,strcmp(datatypes,'DN')} = DN;
    GLMSPopData{n,strcmp(datatypes,'GLMP')} = GLMP;
    
end

save([library2 'GLMSPopData'],'GLMSPopData')
        
% Save variables
set(RFFig.RFpanel,'userdata',rfpanel);
set(RFFig.conpanel,'userdata',conpanel);
set(RFFig.tuningpanel,'userdata',tpanel);
set(gcf,'userdata',RFFig);

end