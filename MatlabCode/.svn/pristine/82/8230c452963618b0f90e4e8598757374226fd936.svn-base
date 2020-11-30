% Mean LinearFilter from White Noise Stimulation of all L and M cones at
% different light levels in pA/R* (equivalent to estimated single photon response)
load('lnmodel_ModelLF.mat')
samplingRate=2e4;
TimeAxis=1:length(Filter);
TimeAxis=TimeAxis./samplingRate;
%%

Guess=[0.6745    0.0216    0.0299    0.5311   34.1814];
% Guess=[2    0.03    0.020  0.3 -58]; %Schanpf and Baylor values

GuessFit=ConeEmpiricalDimFlash(Guess,TimeAxis); 

figure(1)
subplot(2,1,1)
plot(TimeAxis,Filter,'k.')
hold all
plot(TimeAxis,GuessFit,'r-','LineWidth',2)
hold off

% Despite of ggod guess can't get nlinfit to converge. It seems fitting surface is
% pretty flat around here.
DFCoeffs=nlinfit(TimeAxis,Filter,@ConeEmpiricalDimFlash,Guess);
disp(DFCoeffs)
DFModel=ConeEmpiricalDimFlash(DFCoeffs,TimeAxis);
figure(1)
subplot(2,1,2)
plot(TimeAxis,Filter,'k.')
hold all
plot(TimeAxis,DFModel,'r-','LineWidth',2)
hold off