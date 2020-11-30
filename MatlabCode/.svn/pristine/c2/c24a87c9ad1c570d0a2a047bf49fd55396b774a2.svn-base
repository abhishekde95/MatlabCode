% This code is an analysis for GLMP_Hybrid datasets

clear all
close all


% Datasets
%datafile = 'S102111002.nex'; % Symmetric Luminance (Elongated Cartesian grid, non-random sampling)
%datafile = 'S102811002.nex'; % Symmetric Luminance (Elongated Cartesian grid, non-random sampling)
%datafile = 'S103111003.nex'; % Asymmetric Luminance (Elongated Cartesian grid, non-random sampling)
%datafile = 'S103111008.nex'; % Asymmetric Luminance (Elongated Cartesian grid, non-random sampling)
%datafile = 'S110211007.nex'; % Symmetric Luminance (Elongated Cartesian grid, non-random sampling, most data)
%datafile = 'S110411012.nex'; % Asymmetric Col-Lum (Elongated Cartesian grid)  
%datafile = 'S110911003.nex'; % Symmetric Luminance (Elongated Cartesian grid)
%datafile = 'S110911007.nex'; % Symmetric Luminance (Elongated Cartesian grid)
%datafile = 'S111011003.nex'; % Asymmetric Luminance (Elongated Cartesian grid, second to most data)
%datafile = 'S112311006.nex'; % Asymmetric Luminance (Elongated Cartesian grid)

%datafile = 'S112811003.nex'; % Asymmetric Luminance (Elongated Cartesian grid)
%datafile = 'S113011002.nex'; % Asymmetric Col-Lum (Elongated Cartesian grid)
%datafile = 'S120211008.nex'; % Asymmetric M-Cone (Elongated Cartesian grid, crashes in 2D_AllData analysis)
%datafile = 'S120511004.nex'; % Symmetric Luminance (Elongated Cartesian grid)
%datafile = 'S120511008.nex'; % Symmetric Luminance (Elongated Cartesian grid)
%datafile = 'S120811003.nex'; % Symmetric Col-Lum (Elongated radial grid)
%datafile = 'S120911003.nex'; % Asymmetric Col-Lum (Elongated radial grid)
%datafile = 'S121211003.nex'; % Symmetric Luminance (Elongated radial grid)
%datafile = 'S121211006.nex'; % Symmetric Luminance (Elongated radial grid)
%datafile = 'S121511004.nex'; % Symmetric Luminance (Elongated radial grid)

%datafile = 'S011312002.nex'; % Asymmetric Col-Lum (Elongated radial grid)*
%datafile = 'S020112002.nex'; % Symmetric Luminance (Elongated radial grid)
%datafile = 'S020212007.nex'; % Asymmetric Luminance (Elongated radial grid)
%datafile = 'S020812004.nex'; % Asymmetric Col-Lum (Elongated radial grid)*
%datafile = 'S021012007.nex'; % Symmetric Luminance (Elongated radial grid)
%datafile = 'S021712004.nex'; % Symmetric Col-Lum (Elongated radial grid)*
%datafile = 'S021712006.nex'; % Asymmetric Luminance (Elongated radial grid)
%datafile = 'S030112011.nex'; % Asymmetric Luminance (Elongated radial grid)
%datafile = 'S032312002.nex'; % Symmetric Luminance *Plat 2
%datafile = 'S032912004.nex'; % Symmetric Luminance *Plat 2

%datafile = 'S040412002.nex'; % Symmetric Chromatic
%datafile = 'S041012002.nex'; % Symmetric Col-Lum
%datafile = 'S041912004.nex'; % Symmetric Luminance (Logarithmic spacing, bug: theta Unknown)
%datafile = 'S042012002.nex'; % Asymmetric Col-Lum (Possibly Chromatic cell, Logarithmic spacing, bug: Theta Unknown)
%datafile = 'S101812003.nex'; % Asymmetric Luminance (Hybrid)
%datafile = 'S102312007.nex'; % Symmetric Chromatic (Hybrid)
%datafile = 'S102412003.nex'; % M-Cone Responsive (Hybird)
%datafile = 'S110212005.nex'; % Asymmetric Chromatic (Hybrid)


%datafile = 'S110212005.nex'

% Ellipsoidal
%datafile = 'A060112003.nex';% Ellipsoidal Cell! (Apollo's First Cell!)*
%datafile = 'A062612003.nex';% Ellipsoidal Cell!
%datafile = 'A062612004.nex';% Ellipsoidal (Continuation of 062612003. Not getting all waveforms?)
%datafile = 'A070312003.nex';% Ellipsoidal Cell (Tiny Dataset)
%datafile = 'A070312005.nex';% M-Cone Responsive (Too few datapoints...)
%datafile = 'A071012004.nex';% Ellipsoidal Cell (10 repeats)
%datafile = 'A071212005.nex';% Ellipsoidal Cell (Greater chromatic response) (Check iso)
%datafile = 'A071612002.nex';% Ellipsoidal Cell (Check iso)
datafile = 'A071612005.nex';% Ellipsoidal Cell possibly Chromatic?? (Check iso)



%if isdir('Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles')
try
    library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';
    rawdata = nex2stro([char(library) datafile]);
catch
    try isdir('N:\NexFiles\Patrick\')
        library = 'N:\NexFiles\Patrick\';
        rawdata = nex2stro([char(library) datafile]);
    catch
        rawdata = nex2stro;
    end
end

disp(['Now processing ' datafile])
[trial par plat] = OrganizeRawGLMPData(rawdata);


%% Signal Analyses

plat = GLMP_AOV(plat);

plat = GLMP_Response_Variability(plat);

plat = GLMP_TimingAnalysis(plat,1);

%plat = GLMP_ContrastSensitivity(plat);


%% Surface Analyses

%plat = OneDSplineAnalysis(plat);

plat = OneDFit_AllData(plat);

plat = TwoDFit_AllData(plat);

plat = TwoDFit_AllData2(plat,1);

plat = TwoDFit_TwoAxes(plat);

%%
%plat = GLMP_MakeMovie(plat,1);


%% Cone Weight Analyses

%plat = ConeWeightAnalysis(plat);

%plat = ConeWeightAnalysis2(plat);


