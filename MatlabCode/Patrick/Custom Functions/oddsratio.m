function [orat,ci,chi2,p] = oddsratio(a,b,c,d,alpha)
% oddsratio(a,b,c,d,alpha) computes odds ratio and its CI and performs 2x2 contingency table analysis with chi2.
%   Usage: [orat,ci,chi2,p] = oddsratio(a,b,c,d,alpha)
% a,b,c,d are the number of observations that form a 2-by-2 contingency table
% Consider the following 2x2 table.  
% The odds ratio compares the odds of getting disease given exposure
% with the odds of getting disease with no exposure.  The odds is the probabability
% of getting disease divided by 1 minus this probabability.
% The contingency table is in the form
%             Disease
% ex         yes   no
% po    yes   a    b 
% su    no    c    d
% re
%  
% The routine returns the odds ratio (ORAT) and its confidence interval CI 
% at level 1-ALPHA.  The routine also returns a Yates corrected Chi2 value and 
% the p value of accepting the null hypothesis that the odds ratio is 1.
%
% related routines: fisherExact.m


% test vals
% a = 15, b=5, c=9, d=21, alpha = 0.05
% error('forgot to comment out test data')
%

% marginals
n1 = a+b;   % row 1 
n2 = c+d;   % row 2
m1 = a+c;   % col 1
m2 = b+d;   % col 2
n = a+b+c+d;  % total

orat = a*d/(b*c);   % This is the odds ratio

% Confidence interval (Woolf procedure)
q = norminv(1-alpha/2) * sqrt(1./a + 1./b + 1./c + 1./d);
c1 = log(orat) - q;
c2 = log(orat) + q;
ci = exp([c1 c2]);


% Yates corrected chi2
% I got this formula from Bernard & Rosner 4th edition
chi2 = (n * ((abs(a*d - b*c) - n/2).^2) / (m1 * m2)) / (n1 * n2);
p = 1- chi2cdf(chi2,1);

