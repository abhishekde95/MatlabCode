
function [mosaicGLM, goodind] = glmLoad(cellType, varargin)
% Load the parameters for an RGC mosaic measured in an experiment by the
% Chichilnisky Lab, find the average values of each paramter and store a
% mosaic with the average parameters in an isetbio object.
%
% JRG (c) 2016 isetbio team

%% Parse inputs
p = inputParser;
p.addRequired('cellType');
addParameter(p,'cellIndices',   4,     @isnumeric);
addParameter(p,'goodind',    0,     @isnumeric);

p.parse(cellType,varargin{:});

cellType = p.Results.cellType;
cellIndices = p.Results.cellIndices;
goodind  = p.Results.goodind;

%%
% RDT initialization
rdt = RdtClient('isetbio');
switch ieParamFormat(cellType)
    
    % RPE data set - need to put on RDT
    case 'onparasolrpe'
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_RPE_onPar-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/rpe_dataset');
%         data = rdt.readArtifact('mosaicGLM_RPE_onPar', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case 'offparasolrpe'

        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_RPE_offPar-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/rpe_dataset');
%         data = rdt.readArtifact('mosaicGLM_RPE_offPar', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case 'onmidgetrpe'
 
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_RPE_onMid-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/rpe_dataset');
%         data = rdt.readArtifact('mosaicGLM_RPE_onMid', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case 'offmidgetrpe'
 
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_RPE_offMid-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/rpe_dataset');
%         data = rdt.readArtifact('mosaicGLM_RPE_offMid', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
 
    case {'onsbcrpe','sbcrpe'}
     
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_RPE_onSBC-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/rpe_dataset');
%         data = rdt.readArtifact('mosaicGLM_RPE_onSBC', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        
    case{'onparasolapricot'}
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_apricot_ONParasol-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_ONParasol', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case{'offparasolapricot'}
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_apricot_OFFParasol-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_OFFParasol', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case{'onmidgetapricot'}
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_apricot_ONMidget-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_ONMidget', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case{'offmidgetapricot'}
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_apricot_OFFMidget-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_OFFMidget', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case{'onsbcapricot','sbcapricot'}
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_apricot_sbc-1-mat.mat','mosaicGLM');
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_sbc', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
    case 'offparasol'
%         matFileNames = dir([glmFitPath experimentID '/OFF*.mat']);        
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/mosaicGLM_WN_OFFParasol_2013_08_19_6.mat')
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/mosaicGLM_WN_OFFParasol_2013_08_19_6_fits.mat')
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_WN_OFFParasol_2013_08_19_6_fits','mosaicGLM');
%         rdt.crp('resources/data/rgc');
%         data = rdt.readArtifact('mosaicGLM_WN_OFFParasol_2013_08_19_6_fits', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        
%          load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/mosaicGLM_WN_ONParasol_2013_08_19_6_fits.mat')
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/goodind_2013_08_19_6_OFFParasol.mat')
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\goodind_2013_08_19_6_OFFParasol','goodind');
%         rdt = RdtClient('isetbio');
%         rdt.crp('resources/data/rgc');                              
%         data2 = rdt.readArtifact('goodind_2013_08_19_6_OFFParasol', 'type', 'mat');
%         goodind = data2.goodind;
    otherwise % case 'onparasol'
%         matFileNames = dir([glmFitPath experimentID '/ON*.mat']);
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/mosaicGLM_WN_ONParasol_2013_08_19_6_fits.mat')
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\mosaicGLM_WN_ONParasol_2013_08_19_6_fits','mosaicGLM');
%         rdt.crp('resources/data/rgc');
%         data = rdt.readArtifact('mosaicGLM_WN_ONParasol_2013_08_19_6_fits', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;        
        load('C:\Users\Abhishek\Dropbox\UW\Psychophysics_project\RGC_matfiles\goodind_2013_08_19_6_OONParasol','goodind');
%         rdt = RdtClient('isetbio');
%         rdt.crp('resources/data/rgc');
%         data2 = rdt.readArtifact('goodind_2013_08_19_6_ONParasol', 'type', 'mat');
%         goodind = data2.goodind;
end

goodind = 1:length(mosaicGLM);