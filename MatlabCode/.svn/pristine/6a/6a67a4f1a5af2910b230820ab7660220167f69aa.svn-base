function [fit] = ConeEmpiricalDimFlash(coef, t)
% function [fit] = ConeEmpiricalDimFlash(coef, t)
% Following Schnapf, Baylor, 1990: Damped oscillator with S-shaped onset
% but replaced Gaussian by Exponential
% ScFact = coef(1);
% TauR = coef(2); %Rising Phase Time Constant
% TauD = coef(3); %Damping Time Constant
% TauP=coef(4); %Period
% Phi=coef(5); %Phase
% fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
% Created Apr_2010 Angueyra


ScFact = coef(1);% coef(1);
TauR = coef(2); %Rising Phase Time Constant
TauD = coef(3); %Damping Time Constant
TauP = coef(4); %Period
Phi = coef(5); %Phase

% Original Schnapf, Baylor function
% fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD).^2)).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));

% Modified Rieke Angueyra function
fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
