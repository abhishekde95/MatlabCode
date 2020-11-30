function [mu, sigma, initialguesses, stats] = FitROC(data)
% [mu, sigma, initialguesses,stats] = FitROC(data)
%
% Fit an ROC to a buch of points by maximum likelihood
% GDLH 5/23/18
%
% order of columns in data (for standard ROC):
% [n_FAs, n_hits, n_stimabsent, n_stimpresent]
%
% order of columns in data (for nonstandard ROC):
% [n_FAs(laser), n_FAs(no laser), n_stimabsent(laser), n_stimabsent(no laser)]
%
% For non-standard ROC, "fa" = FAs(laser) and "hit" = FAs(no laser)
ALPHA = 0.05; % for confidence interval around muhat

% First getting an initial guess
prop_fa = data(:,1)./data(:,3);
prop_hit =  data(:,2)./data(:,4);
%plot(prop_fa,prop_hit,'k.')

z_fa = 1-norminv(prop_fa);
if any(z_fa < -4 | z_fa > 4)
    keyboard;
end
z_hit = 1-norminv(prop_hit);
if any(z_hit < -4 | z_hit > 4)
    keyboard;
end

b = regress(z_hit, [z_fa ones(size(data,1),1)]);
mu_hat = b(2); % Is this correct?
sigma_hat = 1/b(1);
% Initial guesses
initialguesses = [mu_hat sigma_hat];

options = optimset('Display','iter');
[params, negllik] = fminsearch(@ROCerr,[mu_hat, sigma_hat],options);
mu = params(1);
sigma = params(2);

% Doing hypothesis test on sigma = 1
[params, negllik_nested] = fminsearch(@ROCerr,mu,options);
stats.deviance = -2*(negllik-negllik_nested);
stats.p = 1-chi2cdf(stats.deviance,1);

% Getting a confidence interval on muhat based on likelihood test
bound_err = @(mu) abs(1-chi2cdf(-2*(negllik-ROCerr([mu,sigma])),1)-ALPHA/2);
[upperbound, negllik_tmp1] = fminsearch(bound_err,mu+.01,options);
[lowerbound, negllik_tmp2] = fminsearch(bound_err,mu-.01,options);
stats.muCI = [lowerbound upperbound];

    % ROCerr is a nested function, defined below
    % Inherits "data" because it's nested. Do not modify "data" in this
    % function.
    function total_negllik = ROCerr(parameters)
        mu0 = 0;
        sigma0 = 1;
        mu_candidate = parameters(1);
        if length(parameters) == 2
            sigma_candidate = parameters(2);
        else
            sigma_candidate = 1;
        end
        % preparing for a lookup table for good initial criterion guesses
        tmpcrits1 = linspace(-6,6,10);
        tmpcrits2 = linspace(-6,6,10)*sigma_candidate+mu_candidate;
        generic_ll = @(crit,x,m,s) x(1).*log(1-normcdf(crit,0,1)) + (x(3)-x(1)).*log(normcdf(crit,0,1))+x(2).*log(1-normcdf(crit,m,s)) + (x(4)-x(2)).*log(normcdf(crit,m,s));
        total_negllik = 0;
        for k = 1:size(data,1)
            initialguess1 = interp1(1-normcdf(tmpcrits1,mu0,sigma0),tmpcrits1,data(k,1)/data(k,3));
            initialguess2 = interp1(1-normcdf(tmpcrits2,mu_candidate,sigma_candidate),tmpcrits2, data(k,2)/data(k,4));
            initialguess = nanmean([initialguess1,initialguess2]);
            
            specific_negll = @(x) -1*generic_ll(x,data(k,:),mu_candidate,sigma_candidate);
            [crit_hat, nll] = fminsearch(specific_negll,initialguess);
            total_negllik = total_negllik+nll;
        end
    end
end