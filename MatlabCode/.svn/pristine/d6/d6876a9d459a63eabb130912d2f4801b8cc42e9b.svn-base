% This code creates an interface for GLMS population analyses.
% Created   November, 2014  JPW
function GLMSDPopAnalysis()
global GLMSDPopData

% Set up some variables
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMSD Data\';
end

% Load previously saved data
load([library 'GLMSDPopData.mat'])

% Load in previously analyzed data
data = GLMSDPopData; % #ok<NODEF>

% Create figure
f = figure(456); clf;
set(f,'Pos',[200 200 1000 500])

% Specify other table parameters
colname = {'Filename','Sub','Eccentricity','Stim Size','Psych c50'};
colformat = {'char','numeric','numeric','numeric','numeric'};
colwidth = {'auto' 'auto' 'auto' 'auto' 'auto'};
Popfig.table = uitable('parent',f,'units','normalized','pos',[.01 .01 .58 .98],...
    'data',data(2:end,:),'columnname',colname,...
    'columnformat',colformat,'columnwidth',colwidth,...
    'cellselectioncallback',@cellselect);

% Set up some population analysis buttons
Popfig.poppanel = uipanel('parent',f,'units','normalized',...
    'pos',[.6 .22 .39 .77]);
Popfig.cellpanel = uipanel('parent',f,'units','normalized',...
    'pos',[.6 .01 .39 .2]);

cellpan.newcells = uibutton('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.6 .1 .19 .8],'style','pushbutton','string',{'Look for';'New Cells'},...
    'fontsize',12,'callback',@LookForNewCells);
cellpan.reanalyzeAll = uibutton('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.8 .1 .19 .8],'style','pushbutton','string',{'Reanalyze';'All Cells'},...
    'fontsize',12,'backgroundcolor',[.9 .1 .1],'callback',@ReanalyzeAll);
cellpan.reanalyzeCell = uibutton('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.4 .1 .19 .8],'style','pushbutton','string',{'Reanalyze';'Cell'},...
    'fontsize',12,'callback',@ReanalyzeCell);
cellpan.datafile = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.05 .55 .3 .35],'style','edit','string',[],'fontsize',12);
cellpan.loadConPan = uicontrol('parent',Popfig.cellpanel,'units','normalized',...
    'pos',[.05 .1 .3 .35],'style','pushbutton','string','Load Control Panel',...
    'fontsize',12,'callback',@LoadConPan);

% Population Analyses
poppan.Neurometric = uibutton('parent',Popfig.poppanel,'units','normalized',...
    'pos',[.7 .05 .25 .15],'string',{'Detection' 'Thresholds'},'fontsize',12,...
    'callback',@GLMSDPopGUI_Neurometric);

% Save figure variables
set(Popfig.cellpanel,'userdata',cellpan)
set(Popfig.poppanel,'userdata',poppan)
set(f,'userdata',Popfig)

end


function ReanalyzeCell(~,~)
global GLMSDPopData

% Unpack datafile
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
Popfig = get(1992,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
datafile = get(cellpan.datafile,'string');
sub = cellpan.selectedsub;
idx = cellpan.selectedidx;
load([library 'GLMSDPopData.mat'])

% Load and Organize Rawdata
rawdata = nex2stro([library char(datafile) '.nex']);
[GLMP,~] = OrganizeRawGLMSData(rawdata);

% Re-perform analyses
GLMSGUI_Surface([],[],GLMP,'startanalysis',sub)
SurfFig = get(60,'userdata');
conpan = get(SurfFig.conpanel,'userdata');
fitspanel = get(SurfFig.fitspanel,'userdata');
pval = str2double(get(conpan.paramvals.pval,'string'));
if pval > .05
    GLMSDPopData{idx,3} = '1D';
    GLMSDPopData{idx,4} = pval;
    GLMSDPopData{idx,5} = str2double(get(conpan.paramvals.oneD.rot,'string'));
    parvals = fitspanel.params.oneD;
    surftype = 'surface 7';
    GLMSDPopData{idx,6} = {surftype; parvals};
    GLMSDPopData{idx,7} = fitspanel.conf.oneD;
    
else
    GLMSDPopData{idx,3} = '2D';
    GLMSDPopData{idx,4} = pval;
    GLMSDPopData{idx,5} = str2double(get(conpan.paramvals.twoD.rot,'string'));
    parvals = fitspanel.params.twoD;
    surftype = 'surface 8';
    GLMSDPopData{idx,6} = {surftype; parvals};
    GLMSDPopData{idx,7} = fitspanel.conf.twoD;
end

% Save and display data
Popfig = get(345,'userdata');
dispdata = cat(2,GLMSDPopData(:,1:5),GLMSDPopData(:,7));
set(Popfig.table,'data',dispdata)
save([library 'GLMSDPopData'],'GLMSDPopData')

end


function ReanalyzeAll(~,~)
global GLMSDPopData

disp('Reanalyzing all data...')
disp('(Dont worry, also making a backup.)')

% Set up some variables
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMSD Data/';
end

% Load pop data
load([library 'GLMSDPopData.mat'])

% Save a backup just in case
save([library 'GLMSDPopDataBackup'],'GLMSDPopData')

% Clear out pop data
[GLMSDPopData{2:end,1}] = deal('replace');
save([library 'GLMSDPopData'],'GLMSDPopData')

% Repopulate Pop Data
LookForNewCells

end


function cellselect(a,b)

% Get user data
Popfig = get(345,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
poppan = get(Popfig.poppanel,'userdata');

if isempty(b.Indices)
    return
end
idx = b.Indices(1);
filename = a.Data(idx,1);
set(cellpan.datafile,'string',filename)
cellpan.selectedidx = idx;
cellpan.selectedsub = a.Data{idx,2};

% Find correct library
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/Pics/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\Pics\';
end

try
    A = imread([library filename{1} 's' num2str(cellpan.selectedsub) '.png']);
    axes(poppan.previewWindow)
    image(A)
    set(gca,'XTick',[],'YTick',[])
catch
    disp('No pic to display')
end

% Save user data
set(Popfig.cellpanel,'userdata',cellpan)
set(gcf,'userdata',Popfig)

end


function LoadConPan(~,~)
global GLMP DN

% Get user data
Popfig = get(gcf,'userdata');
cellpan = get(Popfig.cellpanel,'userdata');
datafile = get(cellpan.datafile,'string');

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Unpack datafile
rawdata = nex2stro([char(library) char(datafile) '.nex']);

% Organize Rawdata
[GLMP,DN] = OrganizeRawGLMSDData(rawdata);

% Analysis GUI (interactive data display and analysis)
GLMSD_AnalysisGUI(GLMP,DN);

end


function LookForNewCells(~,~)

if ismac
    library1 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
    library2 = '/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/';
elseif ispc
    library1 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    library2 = 'C:\Users\jpweller\Dropbox\Patrick\GLMSD Data\';
end

% Load in master list of files and previously analyzed data
load([library2 'GLMSDPopData.mat'])
MasterPopData = GLMSDPopData; %#ok<NODEF>
[~,filelist,b] = xlsread([library1  'GLMS_FileNames.xlsx']);

% Compare the two, check for consistency
% if size(GLMSDPopData,1) > size(filelist,1) 
%     GLMSDPopData = GLMSDPopData(2:size(filelist,1)+1,:);
% end
sublist = cat(1,b{:,2});
filesublist = strcat(filelist,num2str(sublist));
%filedata = cat(1,GLMSDPopData(:,1));
%subdata = num2str(cat(1,GLMSDPopData{:,2}));
%filesubdata = strcat(filedata,subdata);
filemaster = cat(1,MasterPopData(2:end,strcmp(MasterPopData(1,:),'Filename')));
submaster = cat(1,MasterPopData{2:end,strcmp(MasterPopData(1,:),'Subunit')});
filesubmaster = strcat(filemaster,num2str(submaster));
counter = 0;

for n = 1:numel(filesublist)
    if any(strcmp(filesublist(n),filesubmaster))
        
        idx = find(strcmp(filesublist(n),filesubmaster),1);
        GLMSDPopData(n,:) = MasterPopData(idx,:);
        
        Popfig = get(1992,'userdata');
        tempdata = GLMSDPopData(:,1:5);
        set(Popfig.table,'data',tempdata)
        save([library2 'GLMSDPopData'],'GLMSDPopData')
        
    else
        
        disp(['Adding ' char(filelist(n)) ' sub # ' num2str(sublist(n)) ' data to the population results...'])
        
        % Add file to the GLMSDPopData.mat file
        datafile = char(filelist(n));
        GLMSDPopData{n+1,strcmp(MasterPopData(1,:),'Filename')} = datafile;
        GLMSDPopData{n+1,strcmp(MasterPopData(1,:),'Subunit')} = sublist(n);
        save([library2 'GLMSDPopData'],'GLMSDPopData')
        
        % Unpack datafile
        rawdata = nex2stro([char(library1) char(datafile) '.nex']);
        
        % Organize Rawdata
        [GLMP,~] = OrganizeRawGLMSData(rawdata);
        [~,rho] = cart2pol(GLMP.rf_x,GLMP.rf_y);
        GLMSDPopData{n+1,strcmp(MasterPopData(1,:),'Eccentricity')} = rho;
        GLMSDPopData{n+1,strcmp(MasterPopData(1,:),'Stim Area')} = numel(GLMP.subunit{sublist(n)}.gridX(1,:));
        
        params = GLMS_ROC(GLMP,sublist(n));
        
        %cellpan.selectedidx = n;
        
        % Save user data
        %set(Popfig.cellpanel,'userdata',cellpan)
        %set(gcf,'userdata',Popfig)
        
        % Save and display data
        Popfig = get(456,'userdata');
        dispdata = GLMSDPopData;
        set(Popfig.table,'data',dispdata)
        save([library2 'GLMSDPopData'],'GLMSDPopData')
        
    end
end

if counter == 0
    disp('No new cells in GLMS_FileNames.xlsx')
else
    disp([num2str(counter) ' new cell(s) added to GLMSDPopData'])
end

end

