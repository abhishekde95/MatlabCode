function GLMSDPopGUI_Neurometric(~,~)
global GLMSPopData

figure(1992); clf;
set (gcf,'pos',[200 200 800 500],'name','Neurometric Analysis');
NMFig.conpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.01 .01 .34 .38]);
NMFig.tuningpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.01 .4 .34 .59]);
NMFig.plotpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.36 .01 .63 .98]);

% Look at 1D or 2D cells
conpanel.nDpanel = uibuttongroup('parent',NMFig.conpanel,'units','normalized',...
    'pos',[.05 .5 .55 .45],'SelectionChangedFcn',@ndselection);
conpanel.oneDbutton = uicontrol(conpanel.nDpanel,'style','radiobutton',...
    'units','normalized','pos',[.1 .65 .8 .3],'string','1D Cells');
conpanel.twoDbutton = uicontrol(conpanel.nDpanel,'style','radiobutton',...
    'units','normalized','pos',[.1 .35 .8 .3],'string','2D Cells');
conpanel.allDbutton = uicontrol(conpanel.nDpanel,'style','radiobutton',...
    'units','normalized','pos',[.1 .05 .8 .3],'string','All Cells');
conpanel.whichnd = 1;

% Look at all directions or selected directions
conpanel.tunedirpanel = uibuttongroup(NMFig.conpanel,'units','normalized',...
    'pos',[.05 .05 .55 .4],'SelectionChangedFcn',@TuneDirSelection);
conpanel.alldirs = uicontrol(conpanel.tunedirpanel,'style','radiobutton',...
    'units','normalized','pos',[.1 .55 .8 .4],'string','All Directions');
conpanel.comparedirs = uicontrol(conpanel.tunedirpanel,'style','radiobutton',...
    'units','normalized','pos',[.1 .05 .8 .5],'string','Compare Directions');
conpanel.whichanal = 'All Directions';

% Buttons for comparing different tune directions
conpanel.binselector = uibuttongroup(NMFig.conpanel,'units','normalized',...
    'pos',[.65 .05 .3 .9],'SelectionChangedFcn',@binID);
conpanel.bin1 = uicontrol(conpanel.binselector,'style','radiobutton',...
    'units','normalized','pos',[.1 .74 .8 .2],'string','Bin # 1','enable','off');
conpanel.bin2 = uicontrol(conpanel.binselector,'style','radiobutton',...
    'units','normalized','pos',[.1 .52 .8 .2],'string','Bin # 2','enable','off');
conpanel.bin3 = uicontrol(conpanel.binselector,'style','radiobutton',...
    'units','normalized','pos',[.1 .28 .8 .2],'string','Bin # 3','enable','off');
conpanel.bin4 = uicontrol(conpanel.binselector,'style','radiobutton',...
    'units','normalized','pos',[.1 .04 .8 .2],'string','Bin # 4','enable','off');
conpanel.whichbin = 1;

% Set up RF Axes
ppanel.axes = axes('parent',NMFig.plotpanel,'units','normalized','pos',[.1 .1 .8 .8]);

% Set up Tuning axes
tpanel.axes = axes('parent',NMFig.tuningpanel,'unit','normalized',...
    'pos',[.1 .1 .8 .8]);

% Pull out some population data
datanames = GLMSPopData(1,:);
nds = GLMSPopData(2:end,strcmp(datanames,'nD'));
tpanel.oneDL = strcmp(nds,'1D');
tpanel.twoDL = strcmp(nds,'2D');
tpanel.allDL = true(size(nds));
tuning = cat(1,GLMSPopData{2:end,strcmp(datanames,'Tuning')});
tpanel.angs1D = tuning(tpanel.oneDL)./180*pi;
tpanel.angs2D = tuning(tpanel.twoDL)./180*pi;
tpanel.angsallD = tuning./180*pi;
rfthetas = cat(1,GLMSPopData{2:end,strcmp(datanames,'RF Theta')});
tpanel.rfthetas1D = rfthetas(tpanel.oneDL);
tpanel.rfthetas2D = rfthetas(tpanel.twoDL);
tpanel.rfthetasallD = rfthetas;
rfrhos = cat(1,GLMSPopData{2:end,strcmp(datanames,'RF Rho')});
tpanel.rfrhos1D = rfrhos(tpanel.oneDL);
tpanel.rfrhos2D = rfrhos(tpanel.twoDL);
tpanel.rfrhosallD = rfrhos;
rfsize = cat(1,GLMSPopData{2:end,strcmp(datanames,'RF Size')});
tpanel.rfsize1D = rfsize(tpanel.oneDL);
tpanel.rfsize2D = rfsize(tpanel.twoDL);
tpanel.rfsizeallD = rfsize;

% Save User Data
set(NMFig.tuningpanel,'userdata',tpanel);
set(NMFig.plotpanel,'userdata',ppanel);
set(NMFig.conpanel,'userdata',conpanel);
set(gcf,'userdata',NMFig);

PlotAllDirs()

end

function binID(~,b)

% Load figure variables
NMFig = get(gcf,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

string = get(b.NewValue,'string');
conpanel.whichbin = str2double(string(end));

% Save User Data
set(NMFig.conpanel,'userdata',conpanel);
set(gcf,'userdata',NMFig);

end

function TuneDirSelection(~,b)

% Load figure variables
NMFig = get(gcf,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Set flag to 'all' or 'select'
conpanel.whichanal = get(b.NewValue,'string');

% Save figure variables
set(NMFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',NMFig);

resetPlots()

end

function ndselection(~,b)

NMFig = get(gcf,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Set flag to 1D or 2D
if strcmp(get(b.NewValue,'string'),'1D Cells')
    conpanel.whichnd = 1;
elseif strcmp(get(b.NewValue,'string'),'2D Cells')
    conpanel.whichnd = 2;
elseif strcmp(get(b.NewValue,'string'),'All Cells')
    conpanel.whichnd = 0;
end

% Save variables
set(NMFig.conpanel,'userdata',conpanel);
set(gcf,'userdata',NMFig);

resetPlots()

end


function resetPlots(~,~)

% Load figure variables
NMFig = get(gcf,'userdata');
conpanel = get(NMFig.conpanel,'userdata');
tpanel = get(NMFig.tuningpanel,'userdata');
ppanel = get(NMFig.plotpanel,'userdata');

if strcmp(conpanel.whichanal,'All Directions')
    
    % Plot RFs
    PlotAllDirs();
    
    % Turn bin buttons off
    set(conpanel.bin1,'enable','off')
    set(conpanel.bin2,'enable','off')
    set(conpanel.bin3,'enable','off')
    set(conpanel.bin4,'enable','off')

elseif strcmp(conpanel.whichanal,'Compare Directions')
    
    % Plot tuning histogram
    axes(tpanel.axes); cla;
    if conpanel.whichnd == 1
        tpanel.hist = rose(tpanel.angs1D);
        %tpanel.hist = rose2(tpanel.angs1D,[],'w');
    elseif conpanel.whichnd == 2
        tpanel.hist = rose(tpanel.angs2D);
    elseif conpanel.whichnd == 0
        tpanel.hist = rose(tpanel.angsallD);
    else
        disp('Porblem with compare directions...')
        keyboard
    end
    set(tpanel.hist,'ButtonDownFcn',@PlotCompareDirs)
    set(tpanel.axes,'ButtonDownFcn',@PlotCompareDirs)

    % Clear rf panel
    axes(ppanel.axes); cla;

    % Turn on bin buttons
    set(conpanel.bin1,'enable','on')
    set(conpanel.bin2,'enable','on')
    set(conpanel.bin3,'enable','on')
    set(conpanel.bin4,'enable','on')
    for n = 1:4
        ppanel.bin{n}.rfthetas = [];
        ppanel.bin{n}.rfrhos = [];
        ppanel.bin{n}.rfsize = [];
        tpanel.bin{n}.angs = [];
    end
    set(tpanel.hist,'color','k','ButtonDownFcn',@PlotCompareDirs)
else
    keyboard
end

% Save figure variables
set(NMFig.tuningpanel,'userdata',tpanel)
set(NMFig.conpanel,'userdata',conpanel)
set(NMFig.plotpanel,'userdata',ppanel)
set(gcf,'userdata',NMFig);

end

function PlotAllDirs(~,~)

% Load Variables
NMFig = get(gcf,'UserData');
tpanel = get(NMFig.tuningpanel,'UserData');
conpanel = get(NMFig.conpanel,'UserData');
ppanel = get(NMFig.plotpanel,'userdata');

if conpanel.whichnd == 1
    angs = tpanel.angs1D;
    %popL = tpanel.oneDL;
    rfthetas = tpanel.rfthetas1D;
    rfrhos = tpanel.rfrhos1D;
    rfsize = tpanel.rfsize1D;
elseif conpanel.whichnd == 2
    angs = tpanel.angs2D;
%    popL = tpanel.twoDL;
    rfthetas = tpanel.rfthetas2D;
    rfrhos = tpanel.rfrhos2D;
    rfsize = tpanel.rfsize2D;
elseif conpanel.whichnd == 0
    angs = tpanel.angsallD;
%    popL = tpanel.allDL;
    rfthetas = tpanel.rfthetasallD;
    rfrhos = tpanel.rfrhosallD;
    rfsize = tpanel.rfsizeallD;
end

% Plot tuning histogram
axes(tpanel.axes);
tpanel.hist = rose2(angs,[],'b');

% Plot RF locations
axes(ppanel.axes); cla;
rfpanel.allbins.plot = plot(rfrhos,rfsize,'k*');
xlabel('Eccentricity')
ylabel('RF Size (# stixels)')

% Save variables
set(NMFig.plotpanel,'userdata',ppanel);
set(NMFig.conpanel,'userdata',conpanel);
set(NMFig.tuningpanel,'userdata',tpanel);
set(gcf,'userdata',NMFig);

end

function PlotCompareDirs(~,~)

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);
[theta,~] = cart2pol(whichpt(1),whichpt(2));
theta = theta+pi;

% Load Variables
NMFig = get(gcf,'UserData');
tpanel = get(NMFig.tuningpanel,'UserData');
conpanel = get(NMFig.conpanel,'UserData');
ppanel = get(NMFig.plotpanel,'userdata');

% 1D or 2D
if conpanel.whichnd == 1
    angs = tpanel.angs1D;
    rfthetas = tpanel.rfthetas1D;
    rfrhos = tpanel.rfrhos1D;
    rfsize = tpanel.rfsize1D;
elseif conpanel.whichnd == 2
    angs = tpanel.angs2D;
    rfthetas = tpanel.rfthetas2D;
    rfrhos = tpanel.rfrhos2D;
    rfsize = tpanel.rfsize2D;
elseif conpanel.whichnd == 0
    angs = tpanel.angsallD;
    rfthetas = tpanel.rfthetasallD;
    rfrhos = tpanel.rfrhosallD;
    rfsize = tpanel.rfsizeallD;
end

% Set up bin colors
cols = ['r';'g';'b';'m'];

% Find the edges of the selected TuneDir
[tout,~] = rose(angs);
edges = unique(tout);
lowerTheta = edges(find(sign(theta-edges)==1,1,'last'))-pi;
upperTheta = edges(find(sign(theta-edges)==-1,1,'first'))-pi;

% Find and plot selected angles
angsIdx = find(angs >= lowerTheta & angs <= upperTheta);
ppanel.bin{conpanel.whichbin}.rfthetas = rfthetas(angsIdx);
ppanel.bin{conpanel.whichbin}.rfrhos = rfrhos(angsIdx);
ppanel.bin{conpanel.whichbin}.rfsize = rfsize(angsIdx);
tpanel.bin{conpanel.whichbin}.angs = angs(angsIdx);

% Pull out greatest eccentricity (to keep graph consistent dimensions)
maxrho = ceil(max(tpanel.rfrhosallD));
maxsize = ceil(max(tpanel.rfsizeallD));

% Plot eccentricity vs size
axes(ppanel.axes); cla; hold on;
for n = 1:4
    %polar(ppanel.bin{n}.rfthetas,ppanel.bin{n}.rfrhos,[cols(n) '*'])
    plot(ppanel.bin{n}.rfrhos,ppanel.bin{n}.rfsize,[cols(n) '*'])
end
set(ppanel.axes,'xlim',[0 maxrho],'ylim',[-0 maxsize])
xlabel('Eccentricity')
ylabel('RF Size (# stixels)')

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
set(NMFig.plotpanel,'userdata',ppanel);
set(NMFig.conpanel,'userdata',conpanel);
set(NMFig.tuningpanel,'userdata',tpanel);
set(gcf,'userdata',NMFig);

end

