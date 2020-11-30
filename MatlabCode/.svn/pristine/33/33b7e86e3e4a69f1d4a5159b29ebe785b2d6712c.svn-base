function [fit] = GovadorvskiiTemplateA1(beta, wavelength)

% Govadorvskii template for A1-based pigments with alpha nad beta band
% Based on Govadorvskii et al., 2000
% Created Dec_2010  Angueyra

A=69.7;
a=0.88;
B=28.0;
b=0.922;
C=-14.9;
c=1.104;
D=0.674;

alphafit=1./(exp(A.*(a-(beta(1)./wavelength)))+exp(B.*(b-(beta(1)./wavelength)))+exp(C.*(c-(beta(1)./wavelength)))+D);


lambdaMaxBeta=189+0.315*beta(1);
b=-40.5+0.195*beta(1); %bandwidth
betafit=beta(2)*exp(-(((wavelength-lambdaMaxBeta)./b).^2));

fit=alphafit+betafit;
