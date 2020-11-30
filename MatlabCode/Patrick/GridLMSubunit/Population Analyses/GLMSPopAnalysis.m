% This code creates an interface for GLMS population analyses.
% Created   November, 2014  JPW
function GLMSPopAnalysis()
global GLMSPopData

% Create figure
f = figure(345); clf;
set(f,'Pos',[200 200 1000 500])

% Set up some population analysis buttons
Popfig.cellpanel = uipanel('parent',f,'units','normalized',...
    'pos',[.7 .01 .29 .48]);
cellpan.datafile = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.03 .675 .45 .3],'style','edit','string',[],'fontsize',12);
cellpan.loadConPan = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.52 .675 .45 .3],'style','pushbutton','string',{'Load Control Panel'},...
    'fontsize',12,'callback',@LoadConPan);
cellpan.analyzeCell = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.03 .35 .45 .3],'style','pushbutton','string','Reanalyze Cell',...
    'fontsize',12,'backgroundcolor',[.7 .3 .3],'callback',@AnalyzeCell);
cellpan.loadOverview = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.52 .35 .45 .3],'style','pushbutton','string','Show Overview',...
    'fontsize',12,'callback',@LoadOverview);
cellpan.newcells = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.03 .025 .45 .3],'style','pushbutton','string',{'Look for New Cells'},...
    'fontsize',12,'callback',@LookForNewCells);
cellpan.reanalyzeAll = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.52 .025 .45 .3],'style','pushbutton','string',{'Reanalyze All Cells'},...
    'fontsize',12,'backgroundcolor',[.9 .1 .1],'callback',@ReanalyzeAll);


% Population Analyses
Popfig.poppanel = uipanel('parent',f,'units','normalized',...
    'pos',[.7 .51 .29 .48]);
poppan.popsurf = uicontrol('parent',Popfig.poppanel,'units','normalized',...
    'pos',[.03 .52 .45 .2],'style','pushbutton','fontsize',12,...
    'string',{'Pop Tuning'},'callback',{@GLMSPopGUI_Tuning});
poppan.params = uicontrol('parent',Popfig.poppanel,'unit','normalized',...
    'pos',[.03 .28 .45 .2],'style','pushbutton','fontsize',12,...
    'string',{'Pop Params'},'callback',{@GLMSPopGUI_Params});
poppan.rfloc = uicontrol('parent',Popfig.poppanel,'units','normalized',...
    'pos',[.03 .04 .45 .2],'style','pushbutton','fontsize',12,...
    'string','RF Locations','callback',@GLMSPopGUI_RFLocations);
poppan.gabors = uicontrol('parent',Popfig.poppanel,'units','normalized',...
    'pos',[.03 .76 .45 .2],'style','pushbutton','fontsize',12,...
    'string','Gabors','callback',@GLMSPopGUI_Gabors);
poppan.neurometric = uicontrol('parent',Popfig.poppanel,'units','normalized',...
    'pos',[.52 .52 .45 .2],'style','pushbutton','fontsize',12,...
    'string',{'Neurometric Surfs'},'callback',@GLMSPopGUI_Neurometric);
poppan.latency = uicontrol('parent',Popfig.poppanel,'unit','normalized',...
    'pos',[.52 .28 .45 .2],'style','pushbutton','fontsize',12,...
    'string','Latency','callback',@GLMSPopGUI_Latency);
poppan.regress = uicontrol('parent',Popfig.poppanel,'unit','normalized',...
    'pos',[.52 .04 .45 .2],'style','pushbutton','fontsize',12,...
    'string','Regression','callback',@GLMSPopGUI_Regress);

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

% Set up population data structure
% poppan.labels{1} = 'Datafile';
% poppan.labels{2} = 'Subunit';
% poppan.labels{3} = 'nD';
% poppan.labels{4} = 'normLLdiff';
% poppan.labels{5} = 'Tuning';
% poppan.labels{6} = 'Tun Std';
% poppan.labels{7} = 'Surface Parameters';
% poppan.labels{8} = 'RF Params';
% poppan.labels{9} = 'NM Params';
% poppan.labels{10} = 'DN';
% poppan.labels{11} = 'GLMP';
% poppan.labels{12} = 'GLMSD';
% poppan.labels{13} = 'GLMPSpikeStats';
% poppan.labels{14} = 'Confidence Intervals';
% poppan.labels{15} = 'Cardinal Hypothesis';
% poppan.labels{16} = 'DN Filters';
poppan.labels = GLMSPopData(1,:);

% Save user data
set(Popfig.cellpanel,'userdata',cellpan)
set(Popfig.poppanel,'userdata',poppan)
set(f,'userdata',Popfig)

% Display data
DisplayTable()

end

function DisplayTable(~,~)
global GLMSPopData

% Load user data
Popfig = get(345,'userdata');

% Pull out variables for display and format datatypes
datanames = GLMSPopData(1,:);
dataexamp = GLMSPopData(2,:);
dispL = false(size(dataexamp));
charL = false(size(dataexamp));
for n = 1:numel(dataexamp)
    if ~isstruct(dataexamp{n}) && size(dataexamp{n},1) == 1
        dispL(n) = true;
        charL(n) = ischar(dataexamp{n});
    end
end
colname = datanames(dispL);
colformat = cell(size(dispL));
colformat(charL) = {'char'};
colformat(~charL) = {'numeric'};
data = GLMSPopData(2:end,dispL); % #ok<NODEF>
%c = uicontextmenu;
%uimenu(c,'Label','Show Overview','Callback',@LoadOverview);
Popfig.table = uitable('parent',345,'units','normalized','pos',[.01 .01 .68 .98],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,... %'columnwidth',colwidth,...
    'cellselectioncallback',@cellselect);

% Save figure variables
set(345,'userdata',Popfig)

end

function LoadOverview(~,~)
global GLMSPopData DN GLMP

% Load user data
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = cellpan.selectedidx+1;
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

function AnalyzeCell(~,~)
global GLMSPopData GLMP DN

% Load figure variables
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
poppan = get(Popfig.poppanel,'userdata');

% Unpack datafile
datafile = get(cellpan.datafile,'string');
if isempty(datafile)
    disp('Must select datafile before reanalysis can occur...')
    return
end
sub = cellpan.selectedsub;
popidx = cellpan.selectedidx+1;
%GLMP = GLMSPopData{popidx,strcmp(poppan.labels,'GLMP')}; %shortcut for quick analysis
%DN = GLMSPopData{popidx,strcmp(poppan.labels,'DN')}; %shortcut for quick analysis

% Analyses (located in subfunctions)
disp(['Analyzing ' char(datafile) ' sub # ' num2str(sub) '...'])
RecRawdata();
%RecRFParams();
%RecSurfParams();
%RecNeurometric();
%RecLatency();
%RecConfInts()

% Record Filename and subunit last, as this is the indicator for reanalysis.
if iscell(datafile)
    GLMSPopData{popidx,strcmp(poppan.labels,'Datafile')} = datafile{1};
else
    GLMSPopData{popidx,strcmp(poppan.labels,'Datafile')} = datafile;
end
GLMSPopData{popidx,strcmp(poppan.labels,'Subunit')} = sub;


    % Load and  Organize Rawdata
    function RecRawdata()
        disp('Unpacking rawdata and saving in population structure...')
        if ismac
            library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
        elseif ispc
            library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
        end
        rawdata = nex2stro([library char(datafile) '.nex']);
        [GLMP,DN] = OrganizeRawGLMSData(rawdata);
        GLMSPopData{popidx,strcmp(poppan.labels,'DN')} = DN;
        GLMSPopData{popidx,strcmp(poppan.labels,'GLMP')} = GLMP;
        %GLMSGUI_Overview();
    end

    % Receptive Field parameters
    function RecRFParams()
        disp('Recording RF Parameters...')
        [theta,rho] = cart2pol(GLMP.rf_x,GLMP.rf_y);
        rf.rfx = GLMP.rf_x;
        rf.rfy = GLMP.rf_y;
        rf.rftheta = theta;
        rf.rfrho = rho;
        rf.dvaperstix = DN.DVAPerStix(1);
        rf.nstix = numel(GLMP.subunit{sub}.gridX{1});
        rf.sqdva = DN.DVAPerStix(1)*rf.nstix;
        GLMSPopData{popidx,strcmp(poppan.labels,'RF Params')} = rf;
    end

    % Surface Fitting
    function RecSurfParams()
        disp('Fitting data with surfaces and recording fitted parameters...')
        
        % Run Surface GUI
        GLMSGUI_Surface([],[],'AnalyzeReq',sub)
        
        % Get parameter values from figure variables
        SurfFig = get(60,'userdata');
        fitspanel = get(SurfFig.fitspanel,'userdata');
        
        % Organize parameters into population structure and save
        oneD.surftype = fitspanel.surftype;
        twoD.surftype = fitspanel.surftype;
        oneD.parvals = fitspanel.params.oneD;
        twoD.parvals = fitspanel.params.twoD;
        oneD.LL = fitspanel.LL.oneD;
        twoD.LL = fitspanel.LL.twoD;
        oneD.Hessian = fitspanel.Hessian.oneD;
        twoD.Hessian = fitspanel.Hessian.twoD;
        oneD.normLL = fitspanel.normLL.oneDnormLL;
        twoD.normLL = fitspanel.normLL.twoDnormLL;
        oneD.rsq = fitspanel.rsq.oneD;
        twoD.rsq = fitspanel.rsq.twoD;
        surfparams = struct('oneD',oneD,'twoD',twoD);
        surfparams.diffnormLL = twoD.normLL - oneD.normLL;
        surfparams.LLRTpval = fitspanel.LLRTpval;
        surfparams.fitcomps = fitspanel.fitcomps;
        Lcc = GLMP.subunit{sub}.uniqueLcc;
        Mcc = GLMP.subunit{sub}.uniqueMcc;
        pred1d = ComputeNakaRushtonJPW(oneD.parvals,[Lcc Mcc],'conicsection');
        pred2d = ComputeNakaRushtonJPW(twoD.parvals,[Lcc Mcc],'conicsection');
        diff1d2d = diff([pred2d pred1d],[],2);
        surfparams.maxdiffnsp = max(abs(diff1d2d));
        GLMSPopData{popidx,strcmp(poppan.labels,'Surface Parameters')} = surfparams;
        
        % Organize into 1D and 2D cells
        GLMSPopData{popidx,strcmp(poppan.labels,'normLLdiff')} = surfparams.diffnormLL;
        if surfparams.diffnormLL < .1 % Hard coded. Would be great to make this variable...
            GLMSPopData{popidx,strcmp(poppan.labels,'nD')} = '1D';
            nd = 'oneD';
        else
            GLMSPopData{popidx,strcmp(poppan.labels,'nD')} = '2D';
            nd = 'twoD';
        end
        GLMSPopData{popidx,strcmp(poppan.labels,'Tuning')} = surfparams.(nd).parvals(end-1)/pi*180;
        GLMSPopData{popidx,strcmp(poppan.labels,'Tun Std')} = sqrt(1/surfparams.(nd).Hessian(end-1,end-1)/pi*180);
    end

    % Neurometric Analysis
    function RecNeurometric()
        disp('Fitting neurometric function...')
        
        % Load neurosurf params that we fit nm with same surftype
        surfparams = GLMSPopData{popidx,strcmp(poppan.labels,'Surface Parameters')};
        
        % 1D fit
        disp('Fitting 1D neurometric function...')
        [parvals,LL,hessian] = GLMSGUI_AUC(sub,surfparams.oneD.surftype,surfparams.oneD.parvals);
        nmparams.oneD.parvals = parvals;
        nmparams.oneD.LL = LL;
        nmparams.oneD.Hessian = hessian;
        nmparams.oneD.surftype = surfparams.oneD.surftype;
        
        % 2D fit
        disp('Fitting 2D neurometric function...')
        [parvals,LL,hessian] = GLMSGUI_AUC(sub,surfparams.twoD.surftype,surfparams.twoD.parvals);
        nmparams.twoD.parvals = parvals;
        nmparams.twoD.LL = LL;
        nmparams.twoD.Hessian = hessian;
        nmparams.twoD.surftype = surfparams.twoD.surftype;
        
        % Save in pop data
        GLMSPopData{popidx,strcmp(poppan.labels,'NM Params')} = nmparams;
    end

    % Latency Analysis
    function RecLatency()
        disp('Calculating Latency...')
        GLMSGUI_GLMPSpikeStats()
        spikefig = get(7548,'userdata');
        scatfig = spikefig.scatfig;
        GLMSPopData{popidx,strcmp(poppan.labels,'GLMPSpikeStats')} = scatfig;
    end

    function RecConfInts()
       disp('Recording Confidence Intervals...')
       GLMSPopGUI_ConfidenceIntervals()
       confintfig = get(55589,'userdata')
    end

% Save and display data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
save([library 'GLMSPopData'],'GLMSPopData')
DisplayTable()
disp(['Finished Reanalyzing Datafile ' GLMP.datafile ' subunit # ' num2str(sub)])

end

function LookForNewCells(~,~)
global GLMSPopData

disp('Adding new cells to the population...')

% Load figure variables
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
poppan = get(Popfig.poppanel,'userdata');
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Load in master list of files and previously analyzed data
[~,~,xlsdata] = xlsread([library  'GLMS_FileNames.xlsx']);
xlsdatafiles = cat(1,xlsdata(:,1));
xlssubs = cat(1,xlsdata{:,2});
xlsdfsubs = strcat(xlsdatafiles,num2str(xlssubs));

% Compare the two, check for consistency
previnfo = GLMSPopData(1,:); % column labels
prevPopData = GLMSPopData; 
prevdfs = cat(1,prevPopData(2:end,strcmp(previnfo,'Datafile')));
prevsubs = num2str(cat(1,prevPopData{2:end,strcmp(previnfo,'Subunit')}));
prevdfsubs = strcat(prevdfs,prevsubs);

% Run through each filename in the spreadsheet. If it's already in the
% population data, transfer the data to the new data structure.  If its
% not, add it.
matchL = nan(numel(xlsdfsubs),1);
for n = 1:numel(xlsdfsubs)
    popidx = n+1;
    if any(strcmp(xlsdfsubs(n),prevdfsubs))
        
        % If datafile was previously analyzed, put it into new pop structure
        prevpopL = find(strcmp(xlsdfsubs(n),prevdfsubs))+1;
        matchL(n) = 1;
        for anal = 1:numel(poppan.labels)
            if any(strcmp(poppan.labels{anal},previnfo))
                GLMSPopData{popidx,anal} = prevPopData{prevpopL,strcmp(poppan.labels{anal},previnfo)};
            end
        end
    else
        matchL(n) = 0;
    end
end
reanalidx = find(~matchL);
for n = 1:numel(reanalidx)
    
    fprintf('\n************************************************\n')
    disp(['Adding datafile ' num2str(n) ' of ' num2str(numel(reanalidx)) ' to the population results...'])
    fprintf('************************************************\n\n')

    % Index into population
    idx = reanalidx(n);
    
    % Mark current datafile in figure variables
    set(cellpan.datafile,'string',(xlsdatafiles(idx)));
    cellpan.selectedsub = xlssubs(idx);
    cellpan.selectedidx = idx;
    
    % Save figure variables before running analysis functions
    set(Popfig.cellpanel,'userdata',cellpan)
    set(345,'userdata',Popfig);
    
    % Run through analyses
    AnalyzeCell();
    
end

DisplayTable()

if numel(reanalidx) == 0
    disp('No new datafiles in GLMS_FileNames.xlsx')
else
    disp([num2str(numel(reanalidx)) ' datafiles added to population.'])
end

save([library 'GLMSPopData'],'GLMSPopData')

end

function ReanalyzeAll(~,~)
global GLMSPopData

disp('Reanalyzing all data...')
disp('(Dont worry, also making a backup.)')

% Set up some variables
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Save a backup just in case
save([library 'GLMSPopDataBackup'],'GLMSPopData')

% Clear out pop data
dfidx = strcmp(GLMSPopData(1,:),'Datafile');
[GLMSPopData{2:end,dfidx}] = deal('replace');
save([library 'GLMSPopData'],'GLMSPopData')

% Repopulate Pop Data
LookForNewCells()

end

function cellselect(a,b)

% Get user data
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');

% Record name and idx of selected datafile
if isempty(b.Indices)
    return
end
idx = b.Indices(1);
filename = a.Data{idx,1};
set(cellpan.datafile,'string',filename);
cellpan.selectedidx = idx;
cellpan.selectedsub = a.Data{idx,2};

% Save figure variables
set(Popfig.cellpanel,'userdata',cellpan)
set(345,'userdata',Popfig)

end

function LoadConPan(~,~)
global GLMP DN GLMSPopData

% Load figure varaibles
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
datafile = get(cellpan.datafile,'string');
if isempty(datafile)
    disp('Must select datafile before Control Panel can be loaded...')
    return
end

% Get saved data structures
popidx = cellpan.selectedidx + 1;
info = GLMSPopData(1,:);
DN = GLMSPopData{popidx,strcmp(info,'DN')};
GLMP = GLMSPopData{popidx,strcmp(info,'GLMP')};

% Analysis GUI (interactive data display and analysis)
GLMS_AnalysisGUI(GLMP,DN);

end

