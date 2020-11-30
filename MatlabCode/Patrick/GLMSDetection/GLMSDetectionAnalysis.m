% This code unpacks and analyzes GLMSDetection datasets.
% Created   July, 2014  JPW
global GLMSD


neuralfile = 'N061513002'; % -L+M
%neuralfile = 'N062713004'; % bipolar chromatic (shitty GLMS STA)
%neuralfile = 'N112713002'; % +L-M
%neuralfile = 'N102413002'; % +L-M
%neuralfile = 'N110514001'; % bipolar luminance
%neuralfile = 'N011215002'; % Horseshoe cell

% Set up some variables
if ismac
    library = ['/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/' neuralfile '/'];
elseif ispc
    library = ['C:\Users\jpweller\Dropbox\Patrick\GLMSD Data\' neuralfile '\'];
end

[~,filelist,~] = xlsread([library neuralfile '.xlsx']);
tempdata = [];
rawdata = [];

for f = 1:size(filelist,1)
    datafile = filelist(f,:);
    if ~strcmp(datafile{1}(1),'%')        
        tempdata = nex2stro([char(library) char(datafile)]);
        rawdata = strocat(rawdata,tempdata);
    end
     
    % Fill in some missing parameters resulting from strocat
    if f == size(filelist,1)
       rawdata.sum.exptParams.rf_x = tempdata.sum.exptParams.rf_x;
       rawdata.sum.exptParams.rf_y = tempdata.sum.exptParams.rf_y;
       rawdata.sum.exptParams.filename = num2str(neuralfile);
       rawdata.sum.exptParams.subunit = tempdata.sum.exptParams.subunit;
    end
    clear tempdata
end


GLMSD = OrganizeRawGLMSDData(rawdata);

GLMSDAnalysisGUI()

ntrperstim = size(GLMSD.Lcc) / size(GLMSD.uniqueLcc);

disp(['Collected ' num2str(ntrperstim) ' trials on each stimulus.'])


