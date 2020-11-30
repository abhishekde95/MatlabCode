% WNmain.m is a script which sends the 'stro' data structure to WNSubunitPlot for further
% processing
tic;
clear all;
close all;
% A = nex2stro('N021915003.nex'); % This file has just 1 subunit
% A = nex2stro('N022715002.nex'); % This file has 2 subunits
A = nex2stro(findfile('N040415002.nex')); % This has 3 subunits
stro = A; % Using only temporarily here for understanding the code WMSubunitPlot.m
WNSubunitPlot(A);
toc;

