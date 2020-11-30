% A new analysis script for GridLMPlane

% Created   4/11/12     JPW
% Modernizing (updating library, parsing analyses into independent codes)  10/29/12    JPW
% Modernizing (All analyses are not in independed code. This script is
%               almost outdated...)   10/31/12   JPW

clear all
close all

library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';

% Bring in raw data
%datafile = 'S102111002.nex'
%datafile = 'S110211007.nex'; % Symmetric Trench (Most Data)
%datafile = 'S110411012_2.nex'; Very noisy.  Possibly not a good dataset.
%datafile = 'S110911003.nex'; % Symmetric Trench
%datafile = 'S110911007.nex'; % Symmetric Trench
%datafile = 'S111011003.nex'; % Asymmetric Trench (Second to Most Data)
%datafile = 'S112311006.nex'; % Asymmetric Trench
%datafile = 'S113011002.nex'; % Complex Cell
%datafile = 'S120511008.nex'; % Symmetric Trench

% Radial Grid
% (Mostly) Luminance
%datafile = 'S121211003.nex';% Symmetric Trench
%datafile = 'S121211006.nex';% Symmetric L-M Trench
%datafile = 'S121511004.nex';% Symmetric L-M Trench (Small Dataset)

% (Mostly) Chromatic
%datafile = 'S041912004.nex';% Chromatic Cell? (Bug: Theta Unknown)

% Ellipsoidal
%datafile = 'S120811003.nex';% Ellipsoidal Cell? Extremely noisy.
%datafile = 'S120911003.nex';% Ellipsoidal Cell? Extremely noisy.
%datafile = 'S042012002.nex';% Ellipsoidal Cell! (Bug: Theta Unknown)
%datafile = 'A060112003.nex';% Ellispoidal Cell! (Apollo's First Cell!)
%datafile = 'A062612003.nex';% Ellispoidal Cell!
%datafile = 'A062612004.nex';% Ellipsoidal (Continuation of 062612003. Not getting all waveforms?)
%datafile = 'A070312003.nex';% Ellipsoidal Cell (Tiny Dataset)
%datafile = 'A070312005.nex';% Ellipsoidal Cell 
%datafile = 'A071012004.nex';% Ellipsoidal Cell (10 repeats)
%datafile = 'A071212005.nex';% Ellipsoidal Cell
%datafile = 'A071612002.nex';% Ellipsoidal Cell
%datafile = 'A071612005.nex';% Ellipsoidal Cell (Same cell as A071612002)

% Plateaus: S-Cone
%datafile = 'S010212003.nex';% Possible Elipse (mostly unmodulated fr)
%datafile = 'S011312002.nex';% Ellipse!
datafile = 'S020112002.nex';% Symmetric L-M Trench (check iso)
%datafile = 'S020212007.nex';% Asymmetric L-M Trench

% Plateaus: 2 Directions of Motion
%datafile = 'S021012007.nex';% Symmetric L-M Trench
%datafile = 'S032312002.nex';% Luminance Cell
%datafile = 'S032912004.nex';% Symmetric L-M Trench
%datafile = 'S040412002.nex';% Chromatic Cell!!
%datafile = 'S041012002.nex';% Ellipsoidal Cell!!
%datafile = 'S040612002.nex';% Great iso, but not much structure...


% Bugs...
%datafile = 'S032812004.nex']);% Dropped 2 headers - must trim orignal file
%datafile = 'S032812006.nex']);% Grating File!

% Setting up figure parameters
fig = 1;
disp(['Now processing ' datafile])
rawdata = nex2stro([char(library) datafile]);
[trial par plat] = OrganizeRawGLMPData(rawdata);


%% 1-D Spline Fit

plat = OneDSplineAnalysis(plat);

%% 2-D Fits

TwoDFit_AllData(plat);

TwoDFit_TwoAxes(plat)

%% Other Analyses

plat = ConeWeightAnalysis(plat);





