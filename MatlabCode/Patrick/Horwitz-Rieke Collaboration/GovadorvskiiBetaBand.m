function [fit] = GovadorvskiiBetaBand(beta,wavelength)

% S_betaband(lambda)=A_betaband.*exp(-(lambda-lambdaMaxBeta)./bandwidth)^2)
% beta=A_betaband (proportion beta to alpha band);
% Based on Govadorvskii et al., 2000
% Created Dec_2010  Angueyra

% lambdaMaxBeta and badwidth depend on lambdaMax of the alpha band
lambdaMaxBeta=189+0.315*beta(1);
b=-40.5+0.195*beta(1); %bandwidth
fit=beta(2)*exp(-(((wavelength-lambdaMaxBeta)./b).^2));

