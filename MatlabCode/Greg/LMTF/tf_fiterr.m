function err = tf_fiterr(params, tfs, sens)
% err = tf_fiterr(params, tfs, sens)
%
% Error function for fitting 1D temporal contrast sensitivity functions.
%
% INPUTS
%   params: {xi, zeta, n, delta_n, tau1, kappa)
%   tfs: temporal frequencies
%   sens: sensitivities (1/threshold)
%
% OUTPUT
%   err: Summed (either squared or absolute) log residuals

if length(params) ~= 6
    error('tf_fiterr requires 6 parameters');
end
xi = params(1);
zeta = params(2);
n = params(3);
delta_n = params(4);
logtau1 = params(5);
kappa = params(6);
logtau2 = logtau1+kappa;
tau1 = 10^logtau1;
tau2 = 10^logtau2;

f1 = (1i*2*pi*tau1.*tfs+1).^-n;
f2 = (1i*2*pi*tau2.*tfs+1).^-(n+delta_n);
f = abs(xi*(f1-zeta*f2));
logres = log10(sens) - log10(f);
err = sum(abs(logres));

end