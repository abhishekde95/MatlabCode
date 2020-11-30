%% This function is intended to optimally match a sin wave with a Gaussian

function [sqerr] = sintoGauss(thetaGuess)

sigma = thetaGuess;

x = linspace(-4*pi,4*pi,10000);
y = x>-pi/2 & x<pi/2;

mu = 0;

gaussian = exp(-(x - mu).^2/(2*sigma.^2));

figure(1); clf; hold on; grid on;
plot(gaussian,'b--')
plot(cos(x).*y,'b')

sqerr = sum((diff([gaussian;cos(x).*y],1,1)).^2)

