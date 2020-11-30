% This code is an analysis for GLMP_Hybrid datasets

library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';

%datafile = 'S101812003.nex'; %H
%datafile = 'S102312007.nex';
%datafile = 'S102412003.nex';
datafile = 'S110212005.nex'

rawdata = nex2stro([char(library) datafile]);
disp(['Now processing ' datafile])
[trial par plat] = OrganizeRawGLMPData(rawdata);

%% 
plat = OneDSplineAnalysis(plat);

TwoDFit_AllData(plat)

TwoDFit_TwoAxes(plat)

%%
plat = ConeWeightAnalysis(plat);
