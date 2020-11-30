function GLMSPopGUI_SpikeStats()

% Set up some variables
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Load previously saved data
load([library 'GLMSPopData.mat'])
[~,filelist,subunits] = xlsread([library  'GLMS_FileNames.xlsx']);
subunits = cat(1,subunits{:,2});

% Set up cell array for data
SpikeMat = cell(size(filelist,1)+1,4);
SpikeMat{1,1} = 'BLmean';
SpikeMat{1,2} = 'BLVar';
SpikeMat{1,3} = 'MaxRMean';
SpikeMat{1,4} = 'MinRMean';


for n = 2:size(filelist,1)+1
    
    % Unpack datafile
    datafile = char(filelist(n-1));
    rawdata = nex2stro([char(library) char(datafile) '.nex']);
    sub = subunits(n-1);
    
    % Organize Rawdata
    [GLMP,~] = OrganizeRawGLMSData(rawdata);
    
    % Populate table
    SpikeMat{n,1} = mean(GLMP.subunit{sub}.blnspikes);
    SpikeMat{n,2} = var(GLMP.subunit{sub}.blnspikes);
    SpikeMat{n,3} = max(GLMP.subunit{sub}.meannspikes);
    SpikeMat{n,4} = min(GLMP.subunit{sub}.meannspikes);
         
end

save([library 'SpikeMat'],'SpikeMat')

keyboard

means = cat(1,SpikeMat{2:end,1})
vars = cat(1,SpikeMat{2:end,2})
mins = cat(1,SpikeMat{2:end,4})



