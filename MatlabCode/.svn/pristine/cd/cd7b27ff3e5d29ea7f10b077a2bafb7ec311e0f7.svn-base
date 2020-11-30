function [x1, x2, m, ps2] = evaluate_posterior(hyp, x, y, minTF, maxTF)
% [x1, x2, m, ps2] = evaluate_posterior(hyp, x, y, minTF, maxTF)
%
% Evaluates a Gaussian process on a hardcoded 50 x 50 meshgridded lattice of 
% points in the plane defined by log10(temporal frequency) (in Hz) and
% color direction in the LM plane (in radians). Returns the mean and the
% variance of the Gaussian process at these 2500 points.
%
%    INPUT
%       hyp: hyperparameters in the standard format (a structure with
%       fields "cov", "lik" and "mean" each of which is a scalar or a
%       column vector)
%       x: nx2 matrix each row of which is [log10(TF) theta]
%       y: log10(threshold)
%       minTF: minimum temporal frequency (log transformed) for meshgrid
%       maxTF: maximum temporal frequency (log transformed) for meshgrid
%
%    OUPUT
%       x1: 50 x 50 matrix of log10(temporal frequency)
%       x2: 50 x 50 matrix of theta (from 0 to pi)
%       m: 2500 x 1 matrix of the mean of the Gaussian process evaluated at
%       (x1, x2)
%       ps2: 2500 x 1 matrix of the variance of the Gaussian process evaluated
%       at (x1, x2)
% 
% Commented by GDLH 4/8/16

npoints = 50;
tf_star= log10(logspace(minTF,maxTF,npoints)); % minTF and maxTF are in log10 units, and so is tf_star
theta_star = linspace(0, pi, npoints+1);
[x1,x2] = meshgrid(tf_star, theta_star(1:end-1)); % x1 is log10(TF), x2 is theta
xs = [x1(:) x2(:)];
[m,~,~,ps2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, xs);