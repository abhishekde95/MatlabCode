function K = covPeriodPi(hyp, x, z, i)

% Stationary covariance function for a smooth periodic function, with period pi 
% in 1d (see covPERiso and covPERard for multivariate data):
%
% k(x,z) = sf^2 * exp( -2*sin^2( pi*||x-z|| )/ell^2 )
%
% where the hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf) ]
%
%
% Based on covPeriodic.m by Carl Edward Rasmussen and Hannes Nickisch, 2011-01-05.
% Modified by GDLH 5/16/14. 
% Removes the hyperparameter 'p' which otherwise gives the kernel
% a flexible period. I want to fix the period at pi. GH
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

[n,D] = size(x);
if D>1, error('Covariance is defined for 1d data only.'), end
ell = exp(hyp(1));
sf2 = exp(2*hyp(2));

% precompute distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sqrt(sq_dist(x'));
  else                                                   % cross covariances Kxz
    K = sqrt(sq_dist(x',z'));
  end
end

if nargin<4                                                        % covariances
    K = sin(K)/ell; K = K.*K; K =   sf2*exp(-2*K);
else                                                               % derivatives
  if i==1
    K = sin(K)/ell; K = K.*K; K = 4*sf2*exp(-2*K).*K;
  elseif i==2
    K = sin(K)/ell; K = K.*K; K = 2*sf2*exp(-2*K);
  else
    error('Unknown hyperparameter')
  end
end