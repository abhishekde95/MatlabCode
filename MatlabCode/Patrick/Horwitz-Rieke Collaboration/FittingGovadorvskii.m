%% Fitting Baylor's Data with Govadorvskii et al. nomograms (Vis NeuroSci,2000, 17, 509-528)
% Rod lsq fit ignores the beta band. Decided to go with parameters derived from fitting in 2 steps, which relies more on falling phase.

% Loading Data
load('BaylorSpectralSens.mat')
Baylor.LCone.Wavelength=BaylorSpectralSens.LCone(:,1)';
Baylor.LCone.SpectralSens=BaylorSpectralSens.LCone(:,2)';
Baylor.LCone.SpectralSens=10.^Baylor.LCone.SpectralSens;
Baylor.LCone.StDev=BaylorSpectralSens.LCone(:,3)';

Baylor.MCone.Wavelength=BaylorSpectralSens.MCone(:,1)';
Baylor.MCone.SpectralSens=BaylorSpectralSens.MCone(:,2)';
Baylor.MCone.SpectralSens=10.^Baylor.MCone.SpectralSens;
Baylor.MCone.StDev=BaylorSpectralSens.MCone(:,3)';

Baylor.SCone.Wavelength=BaylorSpectralSens.SCone(:,1)';
Baylor.SCone.SpectralSens=BaylorSpectralSens.SCone(:,2)';
Baylor.SCone.SpectralSens=10.^Baylor.SCone.SpectralSens;
Baylor.SCone.StDev=BaylorSpectralSens.SCone(:,3)';

Baylor.Rod.Wavelength=BaylorSpectralSens.Rod(:,1)';
Baylor.Rod.SpectralSens=BaylorSpectralSens.Rod(:,2)';
Baylor.Rod.SpectralSens=10.^Baylor.Rod.SpectralSens;
Baylor.Rod.StDev=BaylorSpectralSens.Rod(:,3)';
clear BaylorSpectralSens

Wavelength=[300:1:850]; %Axis for final fits

load('SchnapfHumanRodSpectralSens.mat')
Baylor.HRod.Wavelength=Schnapf.HumanRod(:,1)';
Baylor.HRod.SpectralSens=Schnapf.HumanRod(:,2)';
Baylor.HRod.SpectralSens=10.^Baylor.HRod.SpectralSens;
Baylor.HRod.StDev=Schnapf.HumanRod(:,3)';
clear Schnapf

%% LCones
% Fitting to Govadorvskii Nomogram
LCone.objective='GovadorvskiiTemplateA1';
LCone.x0=[560, 0.14]; %lambdaMax & BetaBandScaling guess
LCone.xdata=Baylor.LCone.Wavelength;
LCone.ydata=Baylor.LCone.SpectralSens;
LCone.lb=[300,0]; %lambdaMax lower bound
LCone.ub=[800,1]; %lambdaMax upper bound
LCone.solver='lsqcurvefit';
LCone.options=optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

fit.LCone=lsqcurvefit(LCone);
disp(fit.LCone);

Baylor.LCone.Fit.LambdaMax=fit.LCone(1);
Baylor.LCone.Fit.BetaBandScaling=fit.LCone(2);
Baylor.LCone.Fit.AlphaBand=GovadorvskiiAlphaBandA1(fit.LCone(1),Wavelength);
Baylor.LCone.Fit.BetaBand=GovadorvskiiBetaBand(fit.LCone,Wavelength);
Baylor.LCone.Nomogram=GovadorvskiiTemplateA1(fit.LCone,Wavelength);

%Final Result
figure(100)
errorbar(Baylor.LCone.Wavelength,log10(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
hold on
plot(Wavelength,log10(Baylor.LCone.Fit.AlphaBand),'b--')
plot(Wavelength,log10(Baylor.LCone.Fit.BetaBand),'c--')
plot(Wavelength,log10(Baylor.LCone.Nomogram),'k-')
xlim([370 840])
ylim([-6 0.1])
title('LCone')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

figure(101)
errorbar(Baylor.LCone.Wavelength,(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
hold on
plot(Wavelength,(Baylor.LCone.Fit.AlphaBand),'b--')
plot(Wavelength,(Baylor.LCone.Fit.BetaBand),'c--')
plot(Wavelength,(Baylor.LCone.Nomogram),'k-')
title('LCone')
xlabel('Wavelength (nm)')
ylabel('Spectral Sensitivity')
hold off

%% MCones
% Fitting to Govadorvskii Nomogram
MCone.objective='GovadorvskiiTemplateA1';
MCone.x0=[535, 0.2]; %lambdaMax & BetaBandScaling guess
MCone.xdata=Baylor.MCone.Wavelength;
MCone.ydata=Baylor.MCone.SpectralSens;
MCone.lb=[300,0]; %lambdaMax lower bound
MCone.ub=[800,1]; %lambdaMax upper bound
MCone.solver='lsqcurvefit';
MCone.options=optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

fit.MCone=lsqcurvefit(MCone);
disp(fit.MCone);

Baylor.MCone.Fit.LambdaMax=fit.MCone(1);
Baylor.MCone.Fit.BetaBandScaling=fit.MCone(2);
Baylor.MCone.Fit.AlphaBand=GovadorvskiiAlphaBandA1(fit.MCone(1),Wavelength);
Baylor.MCone.Fit.BetaBand=GovadorvskiiBetaBand(fit.MCone,Wavelength);
Baylor.MCone.Nomogram=GovadorvskiiTemplateA1(fit.MCone,Wavelength);

%Final Result
figure(102)
errorbar(Baylor.MCone.Wavelength,log10(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
hold on
plot(Wavelength,log10(Baylor.MCone.Fit.AlphaBand),'b--')
plot(Wavelength,log10(Baylor.MCone.Fit.BetaBand),'c--')
plot(Wavelength,log10(Baylor.MCone.Nomogram),'k-')
xlim([370 840])
ylim([-6 0.1])
title('MCone')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

figure(103)
errorbar(Baylor.MCone.Wavelength,(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
hold on
plot(Wavelength,(Baylor.MCone.Fit.AlphaBand),'b--')
plot(Wavelength,(Baylor.MCone.Fit.BetaBand),'c--')
plot(Wavelength,(Baylor.MCone.Nomogram),'k-')
title('MCone')
xlabel('Wavelength (nm)')
ylabel('Spectral Sensitivity')
hold off

%% SCones
% Fitting to Govadorvskii Nomogram
SCone.objective='GovadorvskiiTemplateA1';
SCone.x0=[435, 0.2]; %lambdaMax & BetaBandScaling guess
SCone.xdata=Baylor.SCone.Wavelength;
SCone.ydata=Baylor.SCone.SpectralSens;
SCone.lb=[300,0]; %lambdaMax lower bound
SCone.ub=[800,1]; %lambdaMax upper bound
SCone.solver='lsqcurvefit';
SCone.options=optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

fit.SCone=lsqcurvefit(SCone);
disp(fit.SCone);

Baylor.SCone.Fit.LambdaMax=fit.SCone(1);
Baylor.SCone.Fit.BetaBandScaling=fit.SCone(2);
Baylor.SCone.Fit.AlphaBand=GovadorvskiiAlphaBandA1(fit.SCone(1),Wavelength);
Baylor.SCone.Fit.BetaBand=GovadorvskiiBetaBand(fit.SCone,Wavelength);
Baylor.SCone.Nomogram=GovadorvskiiTemplateA1(fit.SCone,Wavelength);

%Final Result
figure(104)
errorbar(Baylor.SCone.Wavelength,log10(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
plot(Wavelength,log10(Baylor.SCone.Fit.AlphaBand),'g--')
plot(Wavelength,log10(Baylor.SCone.Fit.BetaBand),'c--')
plot(Wavelength,log10(Baylor.SCone.Nomogram),'k-')
xlim([300 650])
ylim([-7 0.1])
title('SCone')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

figure(105)
errorbar(Baylor.SCone.Wavelength,(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
plot(Wavelength,(Baylor.SCone.Fit.AlphaBand),'g--')
plot(Wavelength,(Baylor.SCone.Fit.BetaBand),'c--')
plot(Wavelength,(Baylor.SCone.Nomogram),'k-')
title('SCone')
xlabel('Wavelength (nm)')
ylabel('Spectral Sensitivity')
hold off

%% Rods
% Fitting to Govadorvskii Nomogram
Rod.objective='GovadorvskiiTemplateA1';
Rod.x0=[497.4, 0.3]; %lambdaMax & BetaBandScaling guess
Rod.xdata=Baylor.Rod.Wavelength([2:11]);
Rod.ydata=Baylor.Rod.SpectralSens([2:11]);
Rod.lb=[495,0]; %lambdaMax lower bound
Rod.ub=[503,1]; %lambdaMax upper bound
Rod.solver='lsqcurvefit';
Rod.options=optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

fit.Rod=lsqcurvefit(Rod);
disp(fit.Rod);

Baylor.Rod.Fit.LambdaMax=fit.Rod(1);
Baylor.Rod.Fit.BetaBandScaling=fit.Rod(2);
Baylor.Rod.Fit.AlphaBand=GovadorvskiiAlphaBandA1(fit.Rod(1),Wavelength);
Baylor.Rod.Fit.BetaBand=GovadorvskiiBetaBand(fit.Rod,Wavelength);
Baylor.Rod.Nomogram=GovadorvskiiTemplateA1(fit.Rod,Wavelength);
Baylor.Rod.Nomogram2=GovadorvskiiTemplateA1(Rod.x0,Wavelength);

%Final Result
figure(106)
errorbar(Baylor.Rod.Wavelength,log10(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
hold on
plot(Wavelength,log10(Baylor.Rod.Fit.AlphaBand),'g--')
plot(Wavelength,log10(Baylor.Rod.Fit.BetaBand),'y--')
plot(Wavelength,log10(Baylor.Rod.Nomogram),'b-')
plot(Wavelength,log10(Baylor.Rod.Nomogram2),'c--')
xlim([370 840])
ylim([-6 0.1])
title('Rod')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

figure(107)
errorbar(Baylor.Rod.Wavelength,(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
hold on
plot(Wavelength,(Baylor.Rod.Fit.AlphaBand),'g--')
plot(Wavelength,(Baylor.Rod.Fit.BetaBand),'y--')
plot(Wavelength,(Baylor.Rod.Nomogram),'b-')
plot(Wavelength,(Baylor.Rod.Nomogram2),'c--')
title('Rod')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off
%% Human Rods
% Fitting to Govadorvskii Nomogram
HRod.objective='GovadorvskiiTemplateA1';
HRod.x0=[500, 0.2]; %lambdaMax & BetaBandScaling guess
HRod.xdata=Baylor.HRod.Wavelength;
HRod.ydata=Baylor.HRod.SpectralSens;
HRod.lb=[300,0]; %lambdaMax lower bound
HRod.ub=[800,1]; %lambdaMax upper bound
HRod.solver='lsqcurvefit';
HRod.options=optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

fit.HRod=lsqcurvefit(HRod);
disp(fit.HRod);

Baylor.HRod.Fit.LambdaMax=fit.HRod(1);
Baylor.HRod.Fit.BetaBandScaling=fit.HRod(2);
Baylor.HRod.Fit.AlphaBand=GovadorvskiiAlphaBandA1(fit.HRod(1),Wavelength);
Baylor.HRod.Fit.BetaBand=GovadorvskiiBetaBand(fit.HRod,Wavelength);
Baylor.HRod.Nomogram=GovadorvskiiTemplateA1(fit.HRod,Wavelength);

%Final Result
figure(108)
errorbar(Baylor.HRod.Wavelength,log10(Baylor.HRod.SpectralSens),Baylor.HRod.StDev,'k.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,log10(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'r.','MarkerSize',8)
plot(Wavelength,log10(Baylor.HRod.Fit.AlphaBand),'g--')
plot(Wavelength,log10(Baylor.HRod.Fit.BetaBand),'r--')
plot(Wavelength,log10(Baylor.HRod.Nomogram),'b-')
xlim([370 840])
ylim([-6 0.1])
title('HRod')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

figure(109)
errorbar(Baylor.HRod.Wavelength,(Baylor.HRod.SpectralSens),Baylor.HRod.StDev,'k.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'r.','MarkerSize',8)
plot(Wavelength,(Baylor.HRod.Fit.AlphaBand),'g--')
plot(Wavelength,(Baylor.HRod.Fit.BetaBand),'r--')
plot(Wavelength,(Baylor.HRod.Nomogram),'b-')
title('HRod')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off

%% Summary
% Old Spectra used in Rieke Lab
load('macaque_all.mat')

figure(1000)
subplot(2,2,1)
errorbar(Baylor.SCone.Wavelength,(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
errorbar(Baylor.MCone.Wavelength,(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
errorbar(Baylor.LCone.Wavelength,(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
errorbar(Baylor.HRod.Wavelength,(Baylor.HRod.SpectralSens),Baylor.HRod.StDev,'c.','MarkerSize',8)
plot(Wavelength,(Baylor.SCone.Nomogram),'b-')
plot(Wavelength,(Baylor.Rod.Nomogram),'k-')
plot(Wavelength,(Baylor.MCone.Nomogram),'g-')
plot(Wavelength,(Baylor.LCone.Nomogram),'r-')
plot(Wavelength,(Baylor.HRod.Nomogram),'c-')
title('Fitted Spectra Govadorvskii')
xlabel('Wavelength (nm)')
ylabel('Spectral Sensitivity')
hold off
xlim([300 840])
ylim([-0.2 1.2])

subplot(2,2,2)
errorbar(Baylor.SCone.Wavelength,log10(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,log10(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
errorbar(Baylor.MCone.Wavelength,log10(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
errorbar(Baylor.LCone.Wavelength,log10(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
errorbar(Baylor.HRod.Wavelength,log10(Baylor.HRod.SpectralSens),Baylor.HRod.StDev,'c.','MarkerSize',8)
plot(Wavelength,log10(Baylor.SCone.Nomogram),'b-')
plot(Wavelength,log10(Baylor.Rod.Nomogram),'k-')
plot(Wavelength,log10(Baylor.MCone.Nomogram),'g-')
plot(Wavelength,log10(Baylor.LCone.Nomogram),'r-')
plot(Wavelength,log10(Baylor.HRod.Nomogram),'c-')
title('Fitted Spectra Govadorvskii')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off
xlim([300 840])
ylim([-7 0.1])

subplot(2,2,3)
errorbar(Baylor.SCone.Wavelength,(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
errorbar(Baylor.MCone.Wavelength,(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
errorbar(Baylor.LCone.Wavelength,(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
plot(Spectra(:,1),(Spectra(:,3)),'b-')
plot(Spectra(:,1),(Spectra(:,2)),'k-')
plot(Spectra(:,1),(Spectra(:,4)),'g-')
plot(Spectra(:,1),(Spectra(:,5)),'r-')
title('Fitted Spectra Old')
xlabel('Wavelength (nm)')
ylabel('Spectral Sensitivity')
hold off
xlim([300 840])
ylim([-0.2 1.2])

subplot(2,2,4)
errorbar(Baylor.SCone.Wavelength,log10(Baylor.SCone.SpectralSens),Baylor.SCone.StDev,'b.','MarkerSize',8)
hold on
errorbar(Baylor.Rod.Wavelength,log10(Baylor.Rod.SpectralSens),Baylor.Rod.StDev,'k.','MarkerSize',8)
errorbar(Baylor.MCone.Wavelength,log10(Baylor.MCone.SpectralSens),Baylor.MCone.StDev,'g.','MarkerSize',8)
errorbar(Baylor.LCone.Wavelength,log10(Baylor.LCone.SpectralSens),Baylor.LCone.StDev,'r.','MarkerSize',8)
plot(Spectra(:,1),log10(Spectra(:,3)),'b-')
plot(Spectra(:,1),log10(Spectra(:,2)),'k-')
plot(Spectra(:,1),log10(Spectra(:,4)),'g-')
plot(Spectra(:,1),log10(Spectra(:,5)),'r-')
title('Fitted Spectra Old')
xlabel('Wavelength (nm)')
ylabel('log(Spectral Sensitivity)')
hold off
xlim([300 840])
ylim([-7 0.1])
%% Exporting to Igor
% export=struct;
% export.Wavelength=Wavelength;
% export.SConeSpectralSens=Baylor.SCone.SpectralSens;
% export.MConeSpectralSens=Baylor.MCone.SpectralSens;
% export.LConeSpectralSens=Baylor.LCone.SpectralSens;
% export.RodSpectralSens=Baylor.Rod.SpectralSens;
% export.HumanRodSpectralSens=Baylor.HRod.SpectralSens;
% export.SConeWavelength=Baylor.SCone.Wavelength;
% export.MConeWavelength=Baylor.MCone.Wavelength;
% export.LConeWavelength=Baylor.LCone.Wavelength;
% export.RodWavelength=Baylor.Rod.Wavelength;
% export.HumanRodWavelength=Baylor.HRod.Wavelength;
% export.SConeSD=Baylor.SCone.StDev;
% export.MConeSD=Baylor.MCone.StDev;
% export.LConeSD=Baylor.LCone.StDev;
% export.RodSD=Baylor.Rod.StDev;
% export.HumanRodSD=Baylor.HRod.StDev;
% export.SConeNomogram=Baylor.SCone.Nomogram;
% export.MConeNomogram=Baylor.MCone.Nomogram;
% export.LConeNomogram=Baylor.LCone.Nomogram;
% export.RodNomogram=Baylor.Rod.Nomogram;
% export.HumanRodNomogram=Baylor.HRod.Nomogram;
% export.SConeOldSpectra=Spectra(:,3);
% export.MConeOldSpectra=Spectra(:,4);
% export.LConeOldSpectra=Spectra(:,5);
% export.RodOldSpectra=Spectra(:,2);
% export.OldWavelength=Spectra(:,1);
% 
% options.overwrite = 1;
% eval(sprintf('cd %s',dir.h5_root))
% exportName='FittingGovadorvskii';
% fprintf(1,' - Saving: %s/%s.h5\n',pwd,exportName);
% eval(sprintf('cd %s',dir.root))
% exportStructToHDF5(export,sprintf('%s.h5',exportName),sprintf('%s',exportName),options)
