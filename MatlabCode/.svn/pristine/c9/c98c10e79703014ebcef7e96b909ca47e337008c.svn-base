% This code is for packing up my data to hand to the stats people

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

for n = 67:size(GLMSPopData,1)
    
    disp(n)
    
    % Load and Organize Rawdata
    datafile = GLMSPopData(n+1,1);
    subunit = GLMSPopData{n+1,2};
    rawdata = nex2stro([library char(datafile) '.nex']);
    [GLMP,~] = OrganizeRawGLMSData(rawdata);
    
    % Pull out relevant info
    Lcc = GLMP.subunit{subunit}.Lcc;
    Mcc = GLMP.subunit{subunit}.Mcc;
    nsp = GLMP.subunit{subunit}.nspikes;
    data = cat(2,Lcc,Mcc,nsp);
    
    
    save([library 'N' num2str(n)],'data')
    
    
end
%%

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
CellMat = cell(size(GLMSPopData,1)-1,1);

for n = 1:(size(GLMSPopData,1)-1)
    
    cell = ['N' num2str(n) '.mat'];
    load([library char(cell)])
    CellMat{n} = data;
    
end

save([library 'JPWCellMat'],'CellMat')
